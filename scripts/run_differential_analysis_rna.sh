#!/bin/bash

# How to run
# cd <project_dir>
# bash scripts/run_differential_analysis_rna.sh &> run_differential_analysis_rna.out &
# DO NOT change the name or location of run_differential_analysis_rna.out
# Examples:
# bash scripts/run_differential_analysis_rna.sh &> run_differential_analysis_rna.out &
#
# or to change FDR (padj) of differential analysis:
# bash scripts/run_differential_analysis_rna.sh padj=1 &> run_differential_analysis_rna.out &
#
# or to run and printing all trace commands (i.e. set -x):
# bash scripts/run_create_tracks.sh run=debug &> run_create_tracks.out &

#set -x
set -e

# clear PYTHONPATH so packages are not confused when running the container
unset PYTHONPATH

# get command line arguments
while [[ "$#" -gt 0 ]]; do
	if [[ $1 == "run"* ]];then
		run=$(echo $1 | cut -d '=' -f 2)
		shift
	fi
	if [[ $1 == "padj"* ]];then
		padj=$(echo $1 | cut -d '=' -f 2)
		shift
	fi
	if [[ $1 == "time"* ]];then
		time=$(echo $1 | cut -d '=' -f 2)
		shift
	fi
	if [[ $1 == "batch_adjust"* ]];then
		batch_adjust=$(echo $1 | cut -d '=' -f 2)
		shift
	fi
	if [[ $1 == "help" ]];then
		echo ""
		echo 'usage: bash scripts/run_differential_analysis_rna.sh [OPTION] &> run_differential_analysis_rna.out'
		echo ''
		echo DESCRIPTION
		echo -e '\trun differential RNA-seq analysis'
		echo ''
		echo OPTIONS
		echo ''
		echo help 
		echo -e '\tdisplay this help and exit'
		echo run=echo
		echo -e "\tdo not run, echo all commands. Default is running all commands"
		echo -e "\tif set to "debug", it will run with "set -x""
		echo padj=0.05
		echo -e "\tset FDR of differential genes (as calculated by DESeq2) < 0.05. Default=0.05"
		echo -e "time=1-00:00:00"
		echo -e "\tset SLURM time limit time=DD-HH:MM:SS, where ‘DD’ is days, ‘HH’ is hours, etc."
		echo -e "\tDefault is 1 day."
		echo batch_adjust=yes
		echo -e "\tby default differential analysis was done with replicate batch adjustment."
		echo -e "\tto turn off batch adjustment, set to no.\n"
		exit
	fi
done

date
# converting samples.txt to unix format to remove any invisible extra characters
dos2unix -k samples.txt &> /dev/null

# set default parameters
debug=0
if [[ -z "$run" ]];then
	run=
fi
if [[ -z "$padj" ]];then
	padj=0.05
fi
if [[ -z "$time" ]];then
	time=1-00:00:00
fi
if [[ -z "$batch_adjust" ]];then
        batch_adjust=yes 
fi
if [[ $run == "debug"* ]];then
        set -x
        run=
        debug=1
fi

echo -e "\nRunning differential analysis with $ref_ver as reference. \n"

echo -e "Options used to run:"
echo padj="$padj"
echo time="$time"
echo batch_adjust="$batch_adjust"
echo ""

# project directory
proj_dir=$(pwd)
cd $proj_dir

# this is output files generated in this script will be stored
work_dir=$proj_dir/outputs
echo "all outputs will be stored in $work_dir"
if [ ! -d $work_dir ]; then
	mkdir $work_dir
fi
log_dir=$work_dir/logs
if [ ! -d $log_dir ]; then
	mkdir $log_dir
fi
echo "See $proj_dir/run_differential_analysis_rna.out to check the analysis progress"
echo "All other logs and scripts ran will be stored in $log_dir"
echo ""

# specify number of cpus for featureCounts
ncpus=5
max_cpu=$(lscpu | grep 'CPU(s):' | head -n 1 | awk '{print $2}')
if [[ $max_cpu -lt $ncpus ]]; then
        ncpus=$max_cpu
fi

# singularity image directory
# find based on location of this script
img_dir=$(dirname $(dirname $(readlink -f $0)))

# singularity image name
img_name=rnaseq-pipe-container.sif
# IMPORTANT: It is assumed that:
# scripts to run analysis are in $img_dir/scripts
# reference to run analysis are in $img_dir/ref

echo -e "\nUsing singularity image and scripts in:" ${img_dir} "\n"

# getting SLURM configuration
source $img_dir/scripts/slurm_config_var.sh

# copying scripts ran for records
if [[ ! -d $log_dir/scripts ]];then
	mkdir -p $log_dir/scripts
fi
$(cp $img_dir/scripts/run_differential_analysis_rna.sh $log_dir/scripts)
$(cp $img_dir/scripts/feature_counts_simg.sbatch $log_dir/scripts)
$(cp $img_dir/scripts/run_sartools_simg.sbatch $log_dir/scripts)

# getting samples info from samples.txt
$(sed -e 's/[[:space:]]*$//' samples.txt | sed 's/"*$//g' | sed 's/^"*//g' > samples_tmp.txt)
$(mv samples_tmp.txt samples.txt)
groupname_array=($(awk '!/#/ {print $1}' samples.txt))
repname_array=($(awk '!/#/ {print $3}' samples.txt))
email=$(awk '!/#/ {print $5;exit}' samples.txt | tr -d '[:space:]')
if [[ $email == "NA" ]];then
        email=
fi

# added this to prevent error when using symbolic link that is upstream of
# home directory
# NOTE: for this to work, complete path should be used when creating symbolic link
# example:
# this works:
# cd <project-dir>/ref
# ln -s /gpfs0/home1/gdlessnicklab/lab/data/mm10 .
# this does NOT work:
# # cd <project-dir>/ref
# ln -s /home1/gdlessnicklab/lab/data/mm10 .
# the path should be the path that is returned by 'readlink -f'

# copy run_differential_analysis_rna.out to log_dir
cp $proj_dir/run_differential_analysis_rna.out $log_dir/

#### feature counts ####
# counting number of reads in each feature using subread package featureCounts
# initialize job ids
jid5=
for i in "${!groupname_array[@]}"; do
	groupname=${groupname_array[$i]}
	repname=${repname_array[$i]}
	prefix=${groupname}_${repname}
	cd $proj_dir
	# set -x
	tmp_jid=$(SINGULARITYENV_PYTHONPATH= \
		SINGULARITYENV_run=$run \
		SINGULARITYENV_ncpus=$ncpus \
		SINGULARITYENV_prefix=$prefix \
		SINGULARITYENV_ref_ver=$ref_ver \
		SINGULARITYENV_gtf_file=$gtf_file \
		$run sbatch --output=$log_dir/featureCounts_${prefix}.out \
			--cpus-per-task $ncpus \
			--partition=$high_mem_partition \
			--mail-type=FAIL \
			--mail-user=$email \
			--job-name=featureCounts \
			--time=$time \
			--mem=${med_mem}G \
			--wrap "singularity exec \
				--bind $proj_dir:/mnt \
				--bind $img_dir/scripts:/scripts \
				--bind $gtf_file \
				$img_dir/$img_name \
				/bin/bash /scripts/feature_counts_simg.sbatch"| cut -f 4 -d' ')
	echo "Counting number of reads in each feature for $prefix job id: $tmp_jid"
	echo "See log in $log_dir/featureCounts_${prefix}.out"
	echo ""
	# copy run_differential_analysis_rna.out to log_dir
	cp $proj_dir/run_differential_analysis_rna.out $log_dir/

	if [ -z "$jid5" ]; then
		jid5=$tmp_jid
	else
		jid5=${jid5}:${tmp_jid}
	fi
#set +x
done

cd $proj_dir
# Running DESeq2 analysis via SARTools
# request memory
req_mem=$((${#groupname_array[@]}*50))
if [ "$req_mem" -gt $very_high_mem ]; then
	req_mem=$very_high_mem
fi
req_mem=${req_mem}G
echo "Running differential genes analysis using DESeq2 and SARTools...."
echo "See log is in $log_dir/run_sartools.out"
echo ""
# copy run_differential_analysis_rna.out to log_dir
cp $proj_dir/run_differential_analysis_rna.out $log_dir/

jid6=$(SINGULARITYENV_PYTHONPATH= \
	SINGULARITYENV_run=$run \
	SINGULARITYENV_padj=$padj \
	SINGULARITYENV_email=$email \
	SINGULARITYENV_batch_adjust=$batch_adjust \
	$run sbatch --output=$log_dir/run_sartools.out \
		--job-name=run_sartools \
		--partition=$high_mem_partition \
		--mail-type=FAIL \
		--mail-user=$email \
		--mem=$req_mem \
		--time=$time \
		--dependency=afterok:$jid5 \
		--wrap "singularity exec \
			--bind $proj_dir:/mnt \
			--bind $img_dir/scripts:/scripts \
			$img_dir/$img_name \
				/bin/bash /scripts/run_sartools_simg.sbatch"| cut -f 4 -d' ')

# message if jobs never satisfied or cancelled
check_jid5=$(echo $jid5 | sed 's/:/,/g')
state=($(squeue -j $check_jid5 -h))

while [ ${#state[@]} -ne 0 ];
do
        sleep 10
        state=($(squeue -j $check_jid5 -h))
done

reason=$(squeue -j $jid6 -o "%R" -h)
state=$(sacct -j $jid6 --format=state | tail -n +3 | head -n 1)
if [[ $reason == *"DependencyNeverSatisfied"* || $state == *"CANCELLED"* ]]; then
	scancel $jid6
	echo -e "Feature counts failed. Please check featureCounts_* files in $log_dir\n"
	# copy run_differential_analysis_rna.out to log_dir
	cp $proj_dir/run_differential_analysis_rna.out $log_dir/
	exit
fi

##### re-run final multiqc
jid7=$(SINGULARITYENV_PYTHONPATH= \
	SINGULARITYENV_run=$run \
		SINGULARITYENV_proj_dir=$proj_dir \
		SINGULARITYENV_input_dir=/mnt/outputs \
		$run sbatch --output=$log_dir/multiqc.out \
		--partition=$general_partition \
		--job-name=multiqc \
		--mail-type=FAIL \
		--mail-user=$email \
		--time=$time \
		--dependency=afterok:$jid6 \
		--wrap "singularity exec \
			--bind $proj_dir:/mnt \
			--bind $img_dir/scripts:/scripts \
			$img_dir/$img_name \
			/bin/bash /scripts/multiqc_simg.sbatch"| cut -f 4 -d' ')

# message if jobs never satisfied or cancelled
check_jid6=$(echo $jid6 | sed 's/:/,/g')
state=($(squeue -j $check_jid6 -h))

while [ ${#state[@]} -ne 0 ];
do
        sleep 10
        state=($(squeue -j $check_jid6 -h))
done

reason=$(squeue -j $jid7 -o "%R" -h)
state=$(sacct -j $jid7 --format=state | tail -n +3 | head -n 1)
if [[ $reason == *"DependencyNeverSatisfied"* || $state == *"CANCELLED"* ]]; then
	scancel $jid7
	echo -e "SARTools run failed. Please check run_sartools.out in $log_dir\n"
	# copy run_differential_analysis_rna.out to log_dir
	cp $proj_dir/run_differential_analysis_rna.out $log_dir/
	exit
fi

# delete intermediate bam files
rm -r $proj_dir/outputs/STAR_2pass/Pass1 2> /dev/null || true
rm -r $proj_dir/outputs/STAR_2pass/GenomeForPass2 2> /dev/null || true
rm $proj_dir/outputs/STAR_2pass/Pass2/*Aligned.out.bam 2> /dev/null || true
# delete STAR index only if it is in proj_dir
rm -r $proj_dir/ref/$ref_ver/STAR_index 2> /dev/null || true

message="Differential analysis has been completed\n\
Output files are in $work_dir/diff_analysis_rslt\n\
see $work_dir/diff_analysis_rslt/RNA-seq differential analysis_report.html\n\
for full documentation of differential analysis\n\
see $work_dir/fastqc_rslt for more quality control report\n\n"

jid8=$($run sbatch --dependency=afterok:$jid7 \
		--partition=$general_partition \
		--output=$log_dir/dummy.txt \
		--mail-type=END \
		--mail-user=$email \
		--time=5:00 \
		--job-name=run_differential_analysis_rna \
		--export message="$message",proj_dir=$proj_dir \
		--wrap "echo -e \"$message\"$(date) >> $proj_dir/run_differential_analysis_rna.out"| cut -f 4 -d' ')
# copy run_differential_analysis_rna.out to log_dir
cp $proj_dir/run_differential_analysis_rna.out $log_dir/


#!/bin/bash

# How to run
# cd <project_dir>
# bash /apps/opt/rnaseq-pipeline/scripts/run_differential_analysis_rna.sh &> run_differential_analysis_rna.out &
# Examples:
# bash /apps/opt/rnaseq-pipeline/scripts/run_differential_analysis_rna.sh &> run_differential_analysis_rna.out &
#
# by default, alignment is done to human reference genome hg19 (genome=hg19) unless specified genome=hg38:
# bash /apps/opt/rnaseq-pipeline/scripts/run_differential_analysis_rna.sh genome=hg38 &> run_differential_analysis_rna.out &
#
# or to do nothing but echo all commands:
# bash /apps/opt/rnaseq-pipeline/scripts/run_differential_analysis_rna.sh run=echo &> run_differential_analysis_rna.out &
#
# or to change FDR of differential analysis:
# bash /apps/opt/rnaseq-pipeline/scripts/run_differential_analysis_rna.sh padj=1 \
# &> run_differential_analysis_rna.out &
#
# or to run and printing all trace commands (i.e. set -x):
# bash ~/Steve/virtual_server/cut-n-tag-singularity/scripts/run_create_tracks.sh run=debug &> run_create_tracks.out &

#set -x
set -e
date
# converting samples.txt to unix format to remove any invisible extra characters
dummy=$(dos2unix -k samples.txt)

# clear variable used for optional arguments
unset run time padj
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
	if [[ $1 == "genome"* ]];then
		ref_ver=$(echo $1 | cut -d '=' -f 2)
		shift
	fi
	if [[ $1 == "help" ]];then
		echo 'usage: bash rna-singularity/scripts/run_differential_analysis_rna.sh [OPTION] &> run_differential_analysis_rna.out'
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
		echo -e "if set to "debug", it will run with "set -x""
		echo -e ""
		echo -e "genome=hg19"
		echo -e "\tset reference genome. Default is hg19. Other option: hg38"
		echo -e ""
		echo padj=0.05
		echo -e "\tset FDR of differential genes (as calculated by DESeq2) < 0.05. Default=0.05"
		echo -e "time=<default>"
		echo -e "\tset SLURM time limit time=DD-HH:MM:SS, where ‘DD’ is days, ‘HH’ is hours, etc."
		echo -e "Default=1-00:00:00"
		exit
	fi
done
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
if [[ -z "$ref_ver" ]];then
	ref_ver=hg19
fi
if [[ $run == "debug"* ]];then
        set -x
        run=
        debug=1
fi

# specify version of reference genome options hg19 or hg38
if [[ $ref_ver == 'hg19' ]]; then
	gtf_file=Homo_sapiens.GRCh37.75.gtf
	fasta_file=human_g1k_v37_decoy.fasta
	genome_size=2864785220
elif [[ $ref_ver == 'hg38' ]]; then
	fasta_file=hg38.analysisSet.fa
	gtf_file=hg38.ensGene.gtf
	genome_size=2913022398
fi

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

# specify number of cpus for star, bamCompare and bigwigCompare, featureCounts
ncpus=15

# singularity image directory
# find based on location of this script
img_dir=$(dirname $(dirname $(readlink -f $0)))

# singularity image name
img_name=rnaseq-pipe-container.sif
# IMPORTANT: It is assumed that:
# scripts to run analysis are in $img_dir/scripts
# reference to run analysis are in $img_dir/ref

# copying this script for records
$(cp $img_dir/scripts/run_differential_analysis_rna.sh $log_dir/run_differential_analysis_rna.sh)

# getting samples info from samples.txt
groupname_array=($(awk '!/#/ {print $1}' samples.txt))
repname_array=($(awk '!/#/ {print $3}' samples.txt))
email=$(awk '!/#/ {print $5;exit}' samples.txt | tr -d '[:space:]')
filename_string_array=($(awk '!/#/ {print $6}' samples.txt))
string_pair1_array=($(awk '!/#/ {print $7}' samples.txt))
string_pair2_array=($(awk '!/#/ {print $8}' samples.txt))

#### feature counts ####
# counting number of reads in each feature using subread package featureCounts
# initialize job ids
jid5=
for i in "${!groupname_array[@]}"; do 
	groupname=${groupname_array[$i]}
	repname=${repname_array[$i]}
	string_pair1=${string_pair1_array[$i]}
	string_pair2=${string_pair2_array[$i]}
	filename_string=${filename_string_array[$i]}
	prefix=${groupname}_${repname}
	# 
	cd $proj_dir
	# set -x
	tmp_jid=$(SINGULARITYENV_PYTHONPATH= \
		SINGULARITYENV_run=$run \
		SINGULARITYENV_ncpus=$ncpus \
		SINGULARITYENV_prefix=$prefix \
		SINGULARITYENV_gtf_file=$gtf_file \
		SINGULARITYENV_ref_ver=$ref_ver \
		$run sbatch --output=$log_dir/featureCounts_${prefix}.out \
			--cpus-per-task $ncpus \
			--partition=himem \
			--mail-type=FAIL \
			--mail-user=$email \
			--job-name=feature_ct \
			--time=$time \
			--mem=32G \
			--wrap "singularity exec \
				--bind $proj_dir:/mnt \
				--bind $img_dir/scripts:/scripts \
				--bind $img_dir/ref:/ref \
				$img_dir/$img_name \
				/bin/bash /scripts/feature_counts_simg.sbatch"| cut -f 4 -d' ')
	echo "Counting number of reads in each feature for $prefix job id: $tmp_jid"
	echo "See log in $log_dir/featureCounts_${prefix}.out"
	echo ""
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
if [ "$req_mem" -gt 500 ]; then
	req_mem=500
fi
req_mem=${req_mem}G
echo "Running differential genes analysis using DESeq2 and SARTools...."
echo "See log is in $log_dir/run_sartools.out"
echo ""
jid6=$(SINGULARITYENV_PYTHONPATH= \
	SINGULARITYENV_run=$run \
	SINGULARITYENV_padj=$padj \
	SINGULARITYENV_email=$email \
	$run sbatch --output=$log_dir/run_sartools.out \
		--job-name=run_sartools \
		--partition=himem \
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
	exit
fi	
				
##### re-run final multiqc
jid7=$(SINGULARITYENV_PYTHONPATH= \
	SINGULARITYENV_run=$run \
		SINGULARITYENV_proj_dir=$proj_dir \
		SINGULARITYENV_input_dir=/mnt/outputs \
		$run sbatch --output=$log_dir/multiqc.out \
		--job-name=multiqc \
		--partition=general \
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
	exit
fi	
	
message="Differential analysis has been completed\n\
Output files are in $work_dir/diff_analysis_rslt\n\
see $work_dir/diff_analysis_rslt/RNA-seq differential analysis_report.html\n\
for full documentation of differential analysis\n\
see $work_dir/fastqc_rslt for more quality control report\n\n"

jid8=$($run sbatch --dependency=afterok:$jid7 \
		--output=$log_dir/dummy.txt \
		--mail-type=END \
		--mail-user=$email \
		--time=5:00 \
		--job-name=run_differential_analysis_rna \
		--export message="$message",proj_dir=$proj_dir \
		--wrap "echo -e \"$message\"$(date) >> $proj_dir/run_differential_analysis_rna.out"| cut -f 4 -d' ')

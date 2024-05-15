#!/bin/bash

# print linenumber when set -x is set
export PS4='+${LINENO}:${BASH_SOURCE}: '

# script connecting all the individual scripts to do full rnaseq analysis
# How to run
# cd <project_dir>
# bash scripts/run_rnaseq_full.sh &> run_rnaseq_full.out &
# Examples:
# bash scripts/run_rnaseq_full.sh &> run_rnaseq_full.out &
# to run with specific time limit:
# bash scripts/run_rnaseq_full.sh \
#        	time=DD-HH:MM:SS &> run_rnaseq_full.out &
#
# by default, alignment is done to human reference genome hg19 unless specified using 'genome=hg38':
# bash scripts/run_rnaseq_full.sh genome=hg38 &> run_rnaseq_full.out &
#
# or to do nothing but echo all commands:
# bash scripts/run_rnaseq_full.sh run=echo &> run_rnaseq_full.out &
#
# or to change a bunch of parameters at once:
# bash scripts/run_rnaseq_full.sh run=echo padj=1 time=DD-HH:MM:SS \
# &> run_rnaseq_full.out &
#
# or to run and printing all trace commands (i.e. set -x):
# bash scripts/run_rnaseq_full.sh run=debug &> run_rnaseq_full.out &


#set -x
set -e
# set 'run' to echo to simply echoing all commands
# set to empty to run all commands
# clear variable used for optional arguments
unset run time PYTHONPATH
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
	if [[ $1 == "ref_fa"* ]];then
                ref_fa=$(echo $1 | cut -d '=' -f 2)
                shift
        fi
	if [[ $1 == "ref_gtf"* ]];then
                ref_gtf=$(echo $1 | cut -d '=' -f 2)
                shift
        fi
	if [[ $1 == "batch_adjust"* ]];then
                batch_adjust=$(echo $1 | cut -d '=' -f 2)
                shift
        fi
	if [[ $1 == "ncpus_trim"* ]];then
                ncpus_trim=$(echo $1 | cut -d '=' -f 2)
                shift
        fi
	if [[ $1 == "ncpus_star"* ]];then
                ncpus_star=$(echo $1 | cut -d '=' -f 2)
                shift
        fi
	if [[ $1 == "help" ]];then
		echo ""
		echo 'usage: bash /export/export/apps/opt/rnaseq-pipeline/2.2/scripts/run_rnaseq_full.sh [OPTION] &> run_rnaseq_full.out &'
		echo ''
		echo DESCRIPTION
		echo -e '\trun full RNA-seq analysis: Quality Control, alignment, and differential analysis'
		echo ''
		echo OPTIONS
		echo ''
		echo help 
		echo -e '\tdisplay this help and exit'
		echo run=echo
		echo -e "\tdo not run, echo all commands. Default is running all commands"
		echo -e "\tif set to "debug", it will run with "set -x""
		echo -e "genome=hg19"
		echo -e "\tset reference genome. Default is hg19. Other option: hg38"
		echo -e "\tif using genome other than hg19 or hg38, need to specify both ref_fa and ref_gtf."
                echo -e "ref_fa=/path/to/ref.fa"
                echo -e "\tif using genome other than hg19 or hg38, need to specify ref_fa with path to fasta file"
                echo -e "\tof the reference genome."
                echo -e "ref_gtf=/path/to/ref.gtf"
                echo -e "\tif using genome other than hg19 or hg38, need to specify ref_gtf with path to gtf file"
                echo -e "\tof the reference genome."
		echo padj=0.05
		echo -e "\tset FDR of differential genes (as calculated by DESeq2) < 0.05. Default=0.05"
		echo -e "time=1-00:00:00"
		echo -e "\tset SLURM time limit time=DD-HH:MM:SS, where ‘DD’ is days, ‘HH’ is hours, etc."
		echo -e "\tDefault is 1 day.\n"
		echo batch_adjust=yes
                echo -e "\tby default differential analysis was done with replicate batch adjustment."
                echo -e "\tto turn off batch adjustment, set to no.\n"
		echo ncpus_trim=4
                echo -e "\tset the number of cpus for trimming"
                echo -e "\tDefault is 4 cpus."
                echo -e "\tSet to 1 if pigz error occurs."
		echo ncpus_star=20
                echo -e "\tspecifies the number of cpus used for STAR, bamCompare and feature counts."
                echo -e "\tDefault is 20 cpus."
                echo -e "\tthe number of cpus will be automatically adjusted to max number of cpus if it is"
                echo -e "\tless than the max available cpus."
	exit
	fi
done
date
# converting samples.txt to unix format to remove any invisible extra characters
dos2unix -k samples.txt &> /dev/null

# set default parameters
# note run parameter is set differently here.
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
if [[ -z "$batch_adjust" ]];then
        batch_adjust=yes
fi
if [[ -z "$ncpus_trim" ]];then
        ncpus_trim=4
fi
if [[ -z "$ncpus_star" ]];then
        ncpus_star=20
fi

# singularity image directory
# found based on location of this script
img_dir=$(dirname $(dirname $(readlink -f $0)))

genome_dir=$img_dir/ref/$ref_ver
# set ref_fa to fasta file in genome_dir if variable ref_fa is not defined
if [[ -z $ref_fa ]];then
    	fasta_file=${genome_dir}/$(find $genome_dir -name *.fasta -o -name *.fa | xargs basename)
else
	# else set to ref_fa
        fasta_file=$ref_fa
fi
# set gtf file to ref_gtf in genome_dir if variable ref_gtf is not defined
if [[ -z $ref_gtf ]];then
    	gtf_file=${genome_dir}/$(find $genome_dir -name *.gtf | xargs basename)
else
        gtf_file=$ref_gtf
fi

# run parameter needs to be set differently here
if [[ $run == "debug"* ]];then
        set -x
	run=
	# this is for passing to individual script
	run_debug=debug
fi

# project directory
proj_dir=$(pwd)
cd $proj_dir


echo -e "\nRunning full RNA-seq analysis\n"
echo -e "Options used to run:"
echo padj="$padj"
echo time="$time"
echo genome="$ref_ver"
echo batch_adjust="$batch_adjust"
echo ncpus_trim="$ncpus_trim"
echo ""

echo -e "\nUsing singularity image and scripts in:" ${img_dir} "\n"
# IMPORTANT: It is assumed that:
# scripts to run analysis are in $img_dir/scripts
# reference to run analysis are in $img_dir/ref

work_dir=$proj_dir/outputs
echo -e "All outputs will be stored in $work_dir\n"
if [ ! -d $work_dir ]; then
        mkdir $work_dir
fi
log_dir=$work_dir/logs
if [ ! -d $log_dir ]; then
        mkdir $log_dir
fi
echo -e "Logs and scripts ran will be stored in $log_dir\n"

# copying this script for records
if [[ ! -d $log_dir/scripts ]]; then
	mkdir -p $log_dir/scripts
fi
$(cp $img_dir/scripts/run_rnaseq_full.sh $log_dir/scripts/)

### check and run star index if it doesn't exist
echo ""
echo Check and generate STAR index if genome has not been indexed yet....
echo ""
cp run_rnaseq_full.out $log_dir/

. $img_dir/scripts/run_star_index.sh run=$run_debug time=$time genome=$ref_ver ref_fa=$fasta_file \
	ref_gtf=$gtf_file ncpus_star=$ncpus_star \
	&> run_star_index.out
# also copy to log_dir
cp run_star_index.out $log_dir/
cp run_rnaseq_full.out $log_dir/

# skip checking job if not generating star index.
if [[ $skip_run_star_index == 0 ]];then

	tmp0=$($run sbatch --dependency=$jid0 \
                --time=5:00 \
                --output=$log_dir/dummy_run_star_index.txt \
                --job-name=run_star_index \
                --export message="$message",proj_dir=$proj_dir \
                --wrap "echo -e \"$message\" >> $proj_dir/run_rnaseq_full.out; \
						cp $proj_dir/run_rnaseq_full.out $log_dir/"| cut -f 4 -d' ')
	# message if jobs never satisfied
	check_jid0=$(echo $jid0 | sed 's/:/,/g')
	state=($(squeue -j $check_jid0 -h))

	while [ ${#state[@]} -ne 0 ];
	do
        	sleep 10
        	state=($(squeue -j $check_jid0 -h))
	done

	reason=$(squeue -j $tmp0 -o "%R" -h)
	state=$(sacct -j $tmp0 --format=state | tail -n +3 | head -n 1)
	if [[ $reason == *"DependencyNeverSatisfied"* || $state == *"CANCELLED"* ]]; then
        	scancel $tmp0
        	echo -e "Checking and/or generating STAR index failed. Please check run_star_index.out\n"
        	exit
	fi
fi
echo Done checking and generating STAR index as needed.
echo See run_star_index.out for more details.
echo ""

### trimming and QC
date
echo "Checking whether trimming already ran to completion"
cp $proj_dir/run_rnaseq_full.out $log_dir/
# check if trim_fastqc file exist
# enable nullglob for check_trim
shopt -s nullglob
check_trim=(${proj_dir}/outputs/logs/trim_fastqc_*.out)
# disable nullglob
shopt -u nullglob
if [[ ${#check_trim[@]} -gt 0 ]]; then
	nfastq=$(ls ${proj_dir}/outputs/merged_fastq/*.gz | wc -l)
        # check whether any fail
        ncomplete=$(grep "Analysis complete" ${proj_dir}/outputs/logs/trim_fastqc*.out | wc -l)
        if [[ $ncomplete -eq $nfastq ]];then
        echo "Skip trimming since it's already done."
        else
            	echo Trimming and QC.....see progress in $proj_dir/run_trim_qc.out
                echo ""
		cp $proj_dir/run_rnaseq_full.out $log_dir/
                export run time
                . $img_dir/scripts/run_trim_qc.sh run=$run time=$time ncpus_trim=$ncpus_trim &> run_trim_qc.out
		cp run_trim_qc.out $log_dir/
                echo "Done running trim and QC."
                echo "Read run_trim_qc.out log in ${log_dir} and see whether all steps ran to completion"
                echo ""
		cp $proj_dir/run_rnaseq_full.out $log_dir/
        fi
	cp $proj_dir/run_rnaseq_full.out $log_dir/
else
	echo ""
	echo Trimming and QC.....see progress in run_trim_qc.out
	echo ""
	cp $proj_dir/run_rnaseq_full.out $log_dir/
	cd $proj_dir
	. $img_dir/scripts/run_trim_qc.sh run=$run_debug time=$time ncpus_trim=$ncpus_trim &> run_trim_qc.out
	cp run_trim_qc.out $log_dir/

	message="Done trimming and QC.\n"
	message=${message}"See run_trim_qc.out.\n\n\n"
	message=${message}"Aligning reads to $ref_ver and creating tracks for visualization.....\n"
	message=${message}"See progress in run_align_create_tracks_rna.out\n"

	tmp1=$($run sbatch --dependency=afterok:$jid2 \
		--time=5:00 \
		--output=$log_dir/dummy_run_trim_qc.txt \
		--job-name=run_trim_qc \
		--export message="$message",proj_dir=$proj_dir \
		--wrap "echo -e \"$message\" >> $proj_dir/run_rnaseq_full.out; \
		cp $proj_dir/run_rnaseq_full.out $log_dir/"| cut -f 4 -d' ')
	# message if jobs never satisfied
	check_jid2=$(echo $jid2 | sed 's/:/,/g')
	state=($(squeue -j $check_jid2 -h))

	while [ ${#state[@]} -ne 0 ];
	do
        	sleep 10
        	state=($(squeue -j $check_jid2 -h))
	done

	reason=$(squeue -j $tmp1 -o "%R" -h)
	state=$(sacct -j $tmp1 --format=state | tail -n +3 | head -n 1)
	if [[ $reason == *"DependencyNeverSatisfied"* || $state == *"CANCELLED"* ]]; then
		scancel $tmp1
		echo -e "Trimming and/or QC failed. Please check run_trim_qc.out\n"
	exit
	fi
fi


### aligning reads and creating tracks
cd $proj_dir
#. $img_dir/scripts/run_align_create_tracks_rna.sh run=$run_debug time=$time genome=$ref_ver \
#	ncpus_star=$ncpus_star 2>&1 | tee run_align_create_tracks_rna.out $log_dir/run_align_create_tracks_rna.out
. $img_dir/scripts/run_align_create_tracks_rna.sh run=$run_debug time=$time genome=$ref_ver \
        ncpus_star=$ncpus_star &> $log_dir/run_align_create_tracks_rna.out
cp $log_dir/run_align_create_tracks_rna.out $proj_dir/

message="Done alignment and create tracks for visualization.\n"
message=${message}"See log run_align_create_tracks_rna.out.\n\n\n"
message=${message}"Performing differential genes analysis.....\n"
message=${message}"See progress in run_differential_analysis_rna.out.\n"
cp $proj_dir/run_rnaseq_full.out $log_dir/

tmp=$($run sbatch --dependency=afterok:$jid4c \
		--time=5:00 \
		--output=$log_dir/dummy_run_align_create_tracks_rna.txt \
		--job-name=run_trim_qc \
		--export message="$message",proj_dir=$proj_dir \
		--wrap "echo -e \"$message\" >> $proj_dir/run_rnaseq_full.out; \
		$proj_dir/run_rnaseq_full.out $log_dir/"| cut -f 4 -d' ')
# message if jobs failed
check_jid4c=$(echo $jid4c | sed 's/:/,/g')
state=($(squeue -j $check_jid4c -h))

while [ ${#state[@]} -ne 0 ];
do
        sleep 10
        state=($(squeue -j $check_jid4c -h))
done

reason=$(squeue -j $tmp -o "%R" -h)
state=$(sacct -j $tmp --format=state | tail -n +3 | head -n 1)
if [[ $reason == *"DependencyNeverSatisfied"* || $state == *"CANCELLED"* ]]; then
	scancel $tmp
	echo -e "Alignment and/or track creation failed. Please check run_align_create_tracks_rna.out\n"
	exit
fi
cp $proj_dir/run_rnaseq_full.out $log_dir/

### Running differential genes analysis
cd $proj_dir
. $img_dir/scripts/run_differential_analysis_rna.sh run=$run_debug padj=$padj time=$time genome=$ref_ver batch_adjust=$batch_adjust &> run_differential_analysis_rna.out
cp $proj_dir/run_differential_analysis_rna.out $log_dir/

message="Done differential RNA-seq analysis.\n"
message=$message"See log run_differential_analysis_rna.out\n\n"
message=$message"Done running RNA-seq full analysis\n\n"
message=$message"Output files are in $work_dir\n"
message=$message"Output files for differential analysis are in\n"
message=$message"$work_dir/diff_analysis_rslt\n"
message=$message"RNA-seq_differential_analysis_report.html in \n"
message=$message"$work_dir/diff_analysis_rslt/\n"
message=$message"for full documentation of differential analysis\n"
message=$message"bw_files folder contains the track files for visualization in IGV.\n"
message=$message"STAR_2pass/Pass2/ contains the aligned bam files, and feature counts for each sample\n"
message=$message"fastqc_rslt contains the Quality Control files. See multiqc_report.html for summary of all QC metrics\n"
message=$message"trim folder contains the trimmed fastq files\n\n"
cp $proj_dir/run_rnaseq_full.out $log_dir/

tmp=$($run sbatch --dependency=afterok:$jid8 \
		--time=5:00 \
		--output=$log_dir/dummy_run_differential_analysis_rna.txt \
		--mail-type=END \
		--mail-user=$email \
		--job-name=run_rnaseq_full \
		--export message="$message",proj_dir=$proj_dir \
		--wrap "echo -e \"$message\"$(date) >> $proj_dir/run_rnaseq_full.out; \
		cp $proj_dir/run_rnaseq_full.out $log_dir/"| cut -f 4 -d' ')
# message if jobs failed
check_jid8=$(echo $jid8 | sed 's/:/,/g')
state=($(squeue -j $check_jid8 -h))

while [ ${#state[@]} -ne 0 ];
do
        sleep 10
        state=($(squeue -j $check_jid8 -h))
done

reason=$(squeue -j $tmp -o "%R" -h)
state=$(sacct -j $tmp --format=state | tail -n +3 | head -n 1)
if [[ $reason == *"DependencyNeverSatisfied"* || $state == *"CANCELLED"* ]]; then
	scancel $tmp
	echo -e "Differential RNA-seq analysis failed. Please check run_differential_analysis_rna.out\n"
fi
cp $proj_dir/run_rnaseq_full.out $log_dir/

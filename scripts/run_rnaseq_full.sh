#!/bin/bash

# script connecting all the individual scripts to do full rnaseq analysis
# How to run
# cd <project_dir>
# bash /apps/opt/rnaseq-pipeline/scripts/run_rnaseq_full.sh &> run_rnaseq_full.out &
# Examples:
# bash /apps/opt/rnaseq-pipeline/scripts/run_rnaseq_full.sh &> run_rnaseq_full.out &
# to run with specific time limit:
# 	/apps/opt/rnaseq-pipeline/scripts/run_rnaseq_full.sh \
#        	time=DD-HH:MM:SS &> run_rnaseq_full.out &
# 
# by default, alignment is done to human reference genome hg19 unless specified using 'genome=hg38':
# bash /apps/opt/rnaseq-pipeline/scripts/run_rnaseq_full.sh genome=hg38 &> run_rnaseq_full.out &
# 
# or to do nothing but echo all commands:
# bash /apps/opt/rnaseq-pipeline/scripts/run_rnaseq_full.sh run=echo &> run_rnaseq_full.out &
# 
# or to change a bunch of parameters at once:
# bash /apps/opt/rnaseq-pipeline/scripts/run_rnaseq_full.sh run=echo \
#     padj=1 time=DD-HH:MM:SS \
# &> run_rnaseq_full.out &
# 
# or to run and printing all trace commands (i.e. set -x):
# bash /apps/opt/rnaseq-pipeline/scripts/run_rnaseq_full.sh run=debug &> run_rnaseq_full.out &


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

	if [[ $1 == "help" ]];then
		echo ""
		echo 'usage: bash /apps/opt/rnaseq-pipeline/scripts/run_rnaseq_full.sh [OPTION] &> run_rnaseq_full.out &'
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
		echo -e "\tif using genome other than hg19 or hg38, need to put .fa or .fasta and gtf files"
		echo -e "\tin ref/<genome-name> dir and set genome=<genome-name>."
		echo padj=0.05
		echo -e "\tset FDR of differential genes (as calculated by DESeq2) < 0.05. Default=0.05"
		echo -e "time=1-00:00:00"
		echo -e "\tset SLURM time limit time=DD-HH:MM:SS, where ‘DD’ is days, ‘HH’ is hours, etc."
		echo -e "\tDefault is 1 day.\n"
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
if [[ -z "$ref_ver" ]];then
	ref_ver=hg19
fi
if [[ $run == "debug"* ]];then
        set -x
fi

# project directory
proj_dir=$(pwd)
cd $proj_dir

# singularity image directory
# found based on location of this script
img_dir=$(dirname $(dirname $(readlink -f $0)))

# IMPORTANT: It is assumed that:
# scripts to run analysis are in $img_dir/scripts
# reference to run analysis are in $img_dir/ref

work_dir=$proj_dir/outputs
echo -e "All outputs will be stored in $work_dir\n"
log_dir=$work_dir/logs
echo -e "Logs and scripts ran will be stored in $log_dir\n"

### check and run star index if it doesn't exist
echo ""
echo Check and generate STAR index if genome has not been indexed yet.
echo ""

. $img_dir/scripts/run_star_index.sh run=$run time=$time genome=$ref_ver &> run_star_index.out

# copying this script for records
$(cp $img_dir/scripts/run_rnaseq_full.sh $log_dir/run_rnaseq_full.sh)


### trimming and QC
echo ""
echo Trimming and QC.....see progress in run_trim_qc.out
echo ""

. $img_dir/scripts/run_trim_qc.sh run=$run time=$time &> run_trim_qc.out

message="Done trimming and QC.\n\
See run_trim_qc.out.\n\n\
Aligning reads to $ref_ver and creating tracks for visualization.....\n\
See progress in run_align_create_tracks_rna.out\n\n"

tmp1=$($run sbatch --dependency=afterok:$jid2 \
		--time=5:00 \
		--output=$log_dir/dummy_run_trim_qc.txt \
		--job-name=run_trim_qc \
		--export message="$message",proj_dir=$proj_dir \
		--wrap "echo -e \"$message\" >> $proj_dir/run_rnaseq_full.out"| cut -f 4 -d' ')
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


### aligning reads and creating tracks
cd $proj_dir
. $img_dir/scripts/run_align_create_tracks_rna.sh run=$run time=$time genome=$ref_ver &> run_align_create_tracks_rna.out

message="Done alignment and create tracks for visualization.\n\
See log run_align_create_tracks_rna.out.\n\n\
Performing differential genes analysis .....\n\
See progress in run_differential_analysis_rna.out.\n\n"

tmp=$($run sbatch --dependency=afterok:$jid4c \
		--time=5:00 \
		--output=$log_dir/dummy_run_align_create_tracks_rna.txt \
		--job-name=run_trim_qc \
		--export message="$message",proj_dir=$proj_dir \
		--wrap "echo -e \"$message\" >> $proj_dir/run_rnaseq_full.out"| cut -f 4 -d' ')
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


### Running differential genes analysis
cd $proj_dir
. $img_dir/scripts/run_differential_analysis_rna.sh run=$run padj=$padj time=$time genome=$ref_ver &> run_differential_analysis_rna.out

message="Done differential RNA-seq analysis.\n\
See log run_differential_analysis_rna.out\n\n\
Done running RNA-seq full analysis\n\n\
Output files are in $work_dir\n\
Output files for differential analysis are in\n\
$work_dir/diff_analysis_rslt\n\
see $work_dir/diff_analysis_rslt/RNA-seq \n\
differential analysis_report.html for full documentation of differential analysis\n\
bw_files folder contains the track files for visualization in IGV.\n\
STAR_2pass/Pass2/ contains the aligned bam files, and feature counts for each sample\n\
fastqc_rslt contains the Quality Control files. See multiqc_report.html for summary of all QC metrics\n\
trim folder contains the trimmed fastq files\n\n"

tmp=$($run sbatch --dependency=afterok:$jid8 \
		--time=5:00 \
		--output=$log_dir/dummy_run_differential_analysis_rna.txt \
		--mail-type=END \
		--mail-user=$email \
		--job-name=run_rnaseq_full \
		--export message="$message",proj_dir=$proj_dir \
		--wrap "echo -e \"$message\"$(date) >> $proj_dir/run_rnaseq_full.out"| cut -f 4 -d' ')
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


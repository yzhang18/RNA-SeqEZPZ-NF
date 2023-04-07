#!/bin/bash

# script to trim fastq files and run quality control
# How to run:
# cd <project_dir>
# bash /export/apps/opt/rnaseq-pipeline/2.0/scripts/run_trim_qc.sh &> run_trim_qc.out &
# Examples:
# cd project1
# bash /export/apps/opt/rnaseq-pipeline/2.0/scripts/run_trim_qc.sh &> run_trim_qc.out &
#
# or to run with specific time limit:
# bash /export/apps/opt/rnaseq-pipeline/2.0/scripts/run_trim_qc.sh time=DD-HH:MM:SS  &> run_trim_qc.out &
#
# or to do nothing but echo all commands:
# bash /export/apps/opt/rnaseq-pipeline/2.0/scripts/run_trim_qc.sh run=echo &> run_trim_qc.out &
#
# or to run and printing all trace commands (i.e. set -x):
# bash /export/apps/opt/rnaseq-pipeline/2.0/scripts/run_trim_qc.sh run=debug &> run_trim_qc.out &
#


#set -x
set -e

# clear python path to prevent mixed up of python packages
unset PYTHONPATH
# get command line arguments
while [[ "$#" -gt 0 ]]; do
	if [[ $1 == "run"* ]];then
		run=$(echo $1 | cut -d '=' -f 2)
		shift
	fi
	if [[ $1 == "time"* ]];then
		time=$(echo $1 | cut -d '=' -f 2)
		shift
	fi

	if [[ $1 == "help" ]];then
		echo ''
		echo 'usage: bash /export/apps/opt/rnaseq-pipeline/2.0/scripts/run_trim_qc.sh [OPTION] &> run_trim_qc.out &'
		echo ''
		echo DESCRIPTION
		echo -e '\trun trim and qc for RNA-seq samples'
		echo ''
		echo OPTIONS
		echo ''
		echo help 
		echo -e '\tdisplay this help and exit'
		echo run=echo
		echo -e "\tdo not run, echo all commands. Default is running all commands"
                echo -e "\tif set to "debug", it will run with "set -x""
		echo -e "time=1-00:00:00"
		echo -e "\tset SLURM time limit time=DD-HH:MM:SS, where ‘DD’ is days, ‘HH’ is hours, etc."
		echo -e "\tDefault is 1 day.\n"
		exit
	fi
done

date
echo -e "\nTimming reads for low quality/adapter sequences then run Quality Control\n"
# converting samples.txt to unix format to remove any invisible extra characters
dos2unix -k samples.txt &> /dev/null

# set default parameters
debug=0
if [[ -z "$run" ]];then
	run=
fi
if [[ -z "$time" ]];then
	time=1-00:00:00
fi
if [[ $run == "debug"* ]];then
        set -x
        run=
        debug=1
fi

echo -e "Options used to run:"
echo time="$time"
echo ""

# project directory
proj_dir=$(pwd)
cd $proj_dir

# this is where quality control files generated in this script will be stored
work_dir=$proj_dir/outputs
echo "all outputs will be stored in $work_dir"
if [ ! -d $work_dir ]; then
	mkdir $work_dir
fi
log_dir=$work_dir/logs
if [ ! -d $log_dir ]; then
	mkdir $log_dir
fi
echo "See $proj_dir/run_quality_control.out to check the analysis progress"
echo "All other logs and scripts ran will be stored in $log_dir"
echo ""

# singularity image directory
# find based on location of this script
img_dir=$(dirname $(dirname $(readlink -f $0)))

# singularity image name
img_name=rnaseq-pipe-container.sif
# IMPORTANT: It is assumed that:
# scripts to run analysis are in $img_dir/scripts
# reference to run analysis are in $img_dir/ref

echo -e "\nUsing singularity image and scripts in:" ${img_dir} "\n"

# copying this script for records
$(cp $img_dir/scripts/run_trim_qc.sh $log_dir/run_trim_qc.sh)

# getting samples info from samples.txt
$(sed -e 's/[[:space:]]*$//' samples.txt | sed 's/"*$//g' | sed 's/^"*//g' > samples_tmp.txt)
$(mv samples_tmp.txt samples.txt)
groupname_array=($(awk '!/#/ {print $1}' samples.txt))
repname_array=($(awk '!/#/ {print $3}' samples.txt))
email=$(awk '!/#/ {print $5;exit}' samples.txt | tr -d '[:space:]')
filename_string_array=($(awk '!/#/ {print $6}' samples.txt))
string_pair1_array=($(awk '!/#/ {print $7}' samples.txt))
string_pair2_array=($(awk '!/#/ {print $8}' samples.txt))
# note:
# not sure why the usual echo ${filename_string_array[@]} doesn't work
# but this works
# printf "%s\n" "${filename_string_array[@]}"
# trim adapter and low quality bases

# initialize jobid strings
jid=
tmp_idx=0
cd $proj_dir
# added this to prevent error when using symbolic link that is upstream of
# home directory
# NOTE: for this to work, complete path should be used when creating symbolic link
# example:
# this WORKS:
# cd <project-dir>
# ln -s /gpfs0/home1/gdlessnicklab/lab/data/tmp_data_ct/fastq .
# this does NOT work:
# # cd <project-dir>
# ln -s /home1/gdlessnicklab/lab/data/tmp_data_ct/fastq .
# the path should be the path that is returned by 'readlink -f'

target_link=$(readlink -f fastq)
# note: trim_galore is set to --cores 4 which actually translate to 15 cores
# see trim_galore help for more info
for file in $(find fastq/ -name "*fastq.gz");do
	basefile=$(basename $file)
	idx=$(echo ${filename_string_array[*]} | tr ' ' '\n' | awk -v basefile=$basefile 'basefile ~ $1 {print NR-1}')
	groupname=${groupname_array[$idx]}
	repname=${repname_array[$idx]}
	string_pair1=${string_pair1_array[$idx]}
	string_pair2=${string_pair2_array[$idx]}
	prefix=${groupname}_${repname}

	# need to skip the second read pair 
	if [[ "$basefile" =~ $string_pair1 ]];then
		tmp_idx=$((tmp_idx+1))
		# echo $prefix $basefile $run \
		# $log_dir $email $proj_dir \
		# $groupname $repname $string_pair1 \
		# $string_pair2
		# execute inside singularity
		tmp_jid=$(SINGULARITYENV_PYTHONPATH= \
		SINGULARITYENV_run=$run \
		SINGULARITYENV_prefix=$prefix \
		SINGULARITYENV_file=$file \
		SINGULARITYENV_string_pair1=$string_pair1 \
		SINGULARITYENV_string_pair2=$string_pair2 \
		$run sbatch --output=$log_dir/trim_fastqc_${prefix}_pseudolane${tmp_idx}.out \
			--job-name=trim_fastqc \
			--partition=general \
			--time=$time \
			--mail-type=FAIL \
			--mail-user=$email \
			--cpus-per-task=5 \
			--wrap "singularity exec \
				--bind $img_dir/scripts:/scripts \
				--bind $proj_dir:/mnt \
				--bind $target_link:/fastq \
				$img_dir/$img_name \
				/bin/sh /scripts/trim_fastqc_simg.sbatch" | cut -f 4 -d' ')
		echo "Processing $basefile Job id: $tmp_jid"
		echo "" 
		if [ -z "$jid" ]; then
				jid=$tmp_jid
			else
				jid=${jid}:${tmp_jid}
		fi
		#set +x
	fi
done

##### run multiqc
jid2=$(SINGULARITYENV_PYTHONPATH= \
	SINGULARITYENV_run=$run \
	SINGULARITYENV_proj_dir=$proj_dir \
	SINGULARITYENV_input_dir=/mnt/outputs/fastqc_rslt \
	$run sbatch --output=$log_dir/multiqc.out \
		--job-name=multiqc \
		--partition=general \
		--mail-type=FAIL \
		--mail-user=$email \
		--time=$time \
		--dependency=afterok:$jid \
		--wrap "singularity exec \
			--bind $proj_dir:/mnt \
			--bind $img_dir/scripts:/scripts \
			$img_dir/$img_name \
				/bin/bash /scripts/multiqc_simg.sbatch"| cut -f 4 -d' ')

# message if jobs never satisfied
check_jid=$(echo $jid | sed 's/:/,/g')
state=($(squeue -j $check_jid -h))

while [ ${#state[@]} -ne 0 ];
do
        sleep 10
        state=($(squeue -j $check_jid -h))
done
				
reason=$(squeue -j $jid2 -o "%R" -h)
state=$(sacct -j $jid2 --format=state | tail -n +3 | head -n 1)
if [[ $reason == *"DependencyNeverSatisfied"* || $state == *"CANCELLED"* ]]; then
	scancel $jid2
	echo -e "Either trimming reads or running fastqc failed. Please check trim_fastqc_* files in $log_dir\n"
	exit
else 
	echo ""
	echo "Running multiqc to combine all the quality control files Job id: $jid2"
	echo "$log_dir/multiqc.out contains the commands ran"
	echo ""
fi

message="Done trimming reads and quality control\n\n\
summary of quality control result is in $work_dir/trim/fastqc_rslt/multiqc_report.html\n\
trimmed fastq files are in $work_dir/trim\n\n"

tmp=$($run sbatch --dependency=afterok:$jid2 \
		--output=$log_dir/dummy.txt \
		--time=5:00 \
		--mail-type=END \
		--mail-user=$email \
		--job-name=run_trim_qc \
		--export message="$message",proj_dir=$proj_dir \
		--wrap "echo -e \"$message\"$(date) >> $proj_dir/run_trim_qc.out"| cut -f 4 -d' ')

# message if jobs never satisfied
check_jid2=$(echo $jid2 | sed 's/:/,/g')
state=($(squeue -j $check_jid2 -h))

while [ ${#state[@]} -ne 0 ];
do
        sleep 10
        state=($(squeue -j $check_jid2 -h))
done
				
reason=$(squeue -j $tmp -o "%R" -h)
state=$(sacct -j $tmp --format=state | tail -n +3 | head -n 1)
if [[ $reason == *"DependencyNeverSatisfied"* || $state == *"CANCELLED"* ]]; then
	scancel $tmp
	echo -e "multiqc failed. Please check multiqc.out in $log_dir\n"
fi

#!/bin/bash

# script to create genome STAR index
# How to run:
# cd <my_project_dir>
# bash /export/apps/opt/rnaseq-pipeline/2.0/scripts/run_star_index.sh &> run_star_index.out &
# Examples:
# cd ~/project1
# bash /export/apps/opt/rnaseq-pipeline/2.0/scripts/run_star_index.sh &> run_star_index.out &
#
# or to run with specific time limit:
# bash /export/apps/opt/rnaseq-pipeline/2.0/scripts/run_star_index.sh time=DD-HH:MM:SS  &> run_star_index.out &
#
# by default, alignment is done to human reference genome hg19 unless specified genome=hg38:
# bash /export/apps/opt/rnaseq-pipeline/2.0/scripts/run_star_index.sh genome=hg38 &> run_align.out &
# available genome hg19 or hg38
# if using other genome, genome file (.fasta or fa) and gtf file need to be in <my_project_dir>/ref/<genome.name> dir
# and set as follows: genome=<genome_name> 
#
# or to do nothing but echo all commands:
# bash /export/apps/opt/rnaseq-pipeline/2.0/scripts/run_star_index.sh run=echo &> run_star_index.out &
#
# or to run and printing all trace commands (i.e. set -x):
# bash /export/apps/opt/rnaseq-pipeline/2.0/scripts/run_star_index.sh run=debug &> run_star_index.out &

# set -x
set -e

# clear python path to avoid mix up
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
        if [[ $1 == "genome"* ]];then
                ref_ver=$(echo $1 | cut -d '=' -f 2)
                shift
        fi

        if [[ $1 == "help" ]];then
		echo ""
                echo 'usage:  /export/apps/opt/rnaseq-pipeline/2.0/scripts/run_star_index.sh [OPTION] &> run_star_index.out & '
                echo ''
                echo DESCRIPTION
                echo -e '\tcreating genome index using STAR'
                echo ''
                echo OPTIONS
                echo ''
                echo help
                echo -e '\tdisplay this help and exit'
                echo run=echo
                echo -e "\tdo not run, echo all commands. Default is running all commands"
                echo -e "\tif set to "debug", it will run with "set -x""
                echo -e "genome=hg19"
                echo -e "\tset reference genome. Default is hg19. Other options: hg38"
                echo -e "\tif using genome other than hg19 or hg38, need to put .fa in ref/<genome-name> dir"
                echo -e "\tgenome=<genome-name> and <genome-name>.fa must be in ref/<genome-name>/<genome.name>.fa"
                echo -e "time=1-00:00:00"
                echo -e "\tset SLURM time limit time=DD-HH:MM:SS, where ‘DD’ is days, ‘HH’ is hours, etc."
		echo -e "\tDefault is 1 day."
                echo -e ""
                exit
        fi
done

date
# set default parameters
# variable to turn on tracing i.e. set -x
debug=0
echo $run
if [[ -z "$run" ]];then
        run=
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

# project directory
proj_dir=$(pwd)
cd $proj_dir

# this will be used to print messages as jobs are running
out_file=$proj_dir/run_star_index.out

# specify number of cpus for star, bamCompare and bigwigCompare
ncpus=20
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

echo -e "Generating STAR genome index and get chromosome sizes file.\n"
echo -e "Options used to run:"
echo time="$time"
echo genome="$ref_ver"
echo ""

skip_run_star_index=0
### specify reference genome
# if ref genome is not in $img_dir/ref, set genome_dir to  project dir
genome_dir=$img_dir/ref/$ref_ver
if [[ ! -d $genome_dir ]];then
	genome_dir=$proj_dir/ref/$ref_ver
	if [[ ! -d $genome_dir ]];then
		echo -e "${genome_dir} not found. Please put genome file and gtf file in ${genome_dir}. \n"
	fi
fi
gtf_file=$(find $genome_dir -name *.gtf | xargs basename)
fasta_file=$(find $genome_dir -name *.fasta -o -name *.fa | xargs basename)
# this is where star index will be stored
star_index_dir=$genome_dir/STAR_index
work_dir=$star_index_dir
echo "all outputs will be stored in $work_dir"
# check whether STAR_index exist and not empty
if [ -d "$work_dir" ];then
	if [ -z "$(ls -A $work_dir)" ];then
		echo -e "Generating STAR index.\n"
	else
		# genome index exist, exit script
        	echo -e "run_star_index.sh was not run since genome index already exist.\n"
        	skip_run_star_index=1
	fi
else
	# genome index doesn't exist, run script
	echo -e "Generating STAR index.\n"
fi

log_dir=$proj_dir/outputs/logs
if [ ! -d $log_dir ]; then
	mkdir -p $log_dir
fi

if [[ ! $skip_run_star_index == 1 ]];then
echo "See $proj_dir/run_star_index.out to check the analysis progress"
echo "all other logs will be stored in $log_dir"
echo "log files contain all the commands run"
echo ""

# copying this script for records
$(cp $img_dir/scripts/run_star_index.sh $log_dir/run_star_index.sh)

# converting samples.txt to unix format
dos2unix -k samples.txt &> /dev/null
# getting samples info from samples.txt
# remove trailing tabs, leading and trailing quotations
$(sed -e 's/[[:space:]]*$//' samples.txt | sed 's/"*$//g' | sed 's/^"*//g' > samples_tmp.txt)
$(mv samples_tmp.txt samples.txt)
groupname_array=($(awk '!/#/ {print $1}' samples.txt))
repname_array=($(awk '!/#/ {print $3}' samples.txt))
email=$(awk '!/#/ {print $5;exit}' samples.txt | tr -d '[:space:]')
filename_string_array=($(awk '!/#/ {print $6}' samples.txt))
string_pair1_array=($(awk '!/#/ {print $7}' samples.txt))
string_pair2_array=($(awk '!/#/ {print $8}' samples.txt))

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

# get the "true" path in case it is a symlink
target_link_gtf=$(readlink -f $genome_dir/*.gtf)
target_link_fa=$(readlink -f $genome_dir/*.fa*)
target_fa_name=$(basename $target_link_fa)
target_gtf_name=$(basename $target_link_gtf)
target_fa_dir=$(dirname $target_link_gtf)
target_gtf_dir=$(dirname $target_link_gtf)

#### generate index ####
jid0=$(SINGULARITYENV_PYTHONPATH= \
	SINGULARITYENV_run=$run \
		SINGULARITYENV_ncpus=$ncpus \
		SINGULARITYENV_prefix=$prefix \
		SINGULARITYENV_ref_ver=$ref_ver \
		SINGULARITYENV_target_fa_name=$target_fa_name \
		SINGULARITYENV_target_gtf_name=$target_gtf_name \
		$run sbatch --output=$log_dir/star_index.out \
			--cpus-per-task $ncpus \
			--partition=himem \
			--mail-type=FAIL \
			--mail-user=$email \
			--job-name=star_index \
			--time=$time \
			--wrap "singularity exec \
				--bind $proj_dir:/mnt \
				--bind $img_dir/scripts:/scripts \
				--bind $genome_dir:/ref \
				--bind $target_fa_dir:/ref_fa \
				--bind $target_gtf_dir:/ref_gtf \
				$img_dir/$img_name \
				/bin/bash /scripts/star_index_simg.sbatch"| cut -f 4 -d' ')
echo "Generating STAR index job id: $jid0"
echo ""

# check to make sure jobs are completed. Print messages if not.
msg_ok="run_star_index.sh completed successfully.\n"
msg_ok="${msg_ok}STAR index files are in ${genome_dir}/STAR_index.\n"
msg_fail="One of the steps in run_star_index.sh failed\n"
jid_to_check=$jid0
check_run_star_index_jid=$($run sbatch \
        --output=$log_dir/check_run_star_index.out \
        --mail-type=END \
        --mail-user=$email \
        --wait \
        --time=$time \
        --parsable \
        --job-name=check_run_star_index \
        --export=out_file="$out_file",jid_to_check="$jid_to_check",msg_ok="$msg_ok",msg_fail="$msg_fail",debug="$debug" \
	--wrap "bash $img_dir/scripts/check_job.sh")
fi # skip run_star_index
cp $proj_dir/samples.txt $log_dir/

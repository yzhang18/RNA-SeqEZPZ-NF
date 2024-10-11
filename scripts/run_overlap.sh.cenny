#!/bin/bash

# script to run shiny app
# How to run
# cd <project_dir>
# bash /export/export/apps/opt/rnaseq-pipeline/2.2/scripts/run_overlap.sh
# Examples:
# bash /export/export/apps/opt/rnaseq-pipeline/2.2/scripts/run_overlap.sh
# <project_dir> is where the fastq and outputs directory are
# project directory where the fastq and outputs directory are
#
# to run with all commands printed
# bash /export/export/apps/opt/rnaseq-pipeline/2.2/scripts/run_overlap.sh run=debug
#
# Default time limit is 1 day to change to 1 day, 2 hours, 20 minutes and 30 sec
# bash /export/export/apps/opt/rnaseq-pipeline/2.2/scripts/run_overlap.sh time=1-2:20:30

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
		echo ""
                echo 'usage: bash /export/export/apps/opt/rnaseq-pipeline/2.2/scripts/run_overlap.sh [OPTION]'
                echo ''
                echo DESCRIPTION
                echo -e '\trun shiny app to show overlaps'
                echo ''
                echo OPTIONS
                echo ''
                echo help
                echo -e '\tdisplay this help and exit'
                echo run=echo
                echo -e "\tdo not run, echo all commands. Default is running all commands"
                echo -e "\tif set to "debug", it will run with "set -x""
		echo -e ""
                echo -e "time=1-00:00:00"
                echo -e "\tset SLURM time limit time=DD-HH:MM:SS, where ‘DD’ is days, ‘HH’ is hours, etc."
                echo -e ""
                exit
        fi
done

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

proj_dir=$(pwd)
cd $proj_dir

work_dir=$proj_dir/outputs

# singularity image directory
# find based on location of this script
img_dir=$(dirname $(dirname $(readlink -f $0)))

# singularity image name
img_name=rnaseq-pipe-container.sif
# IMPORTANT: It is assumed that:
# scripts to run analysis are in $img_dir/scripts
# reference to run analysis are in $img_dir/ref

echo -e "\nUsing singularity image and scripts in:" ${img_dir} "\n"

echo -e "Options used to run:"
echo time="$time"
echo ""

# find available port
FROM=1035
TO=65535
port_num=$(singularity exec $img_dir/$img_name comm -23 \
<(seq "$FROM" "$TO" | sort) \
<(ss -Htan | awk '{print $4}' | cut -d':' -f2 | sort -u) \
| shuf | head -n 1)

# Activate environment where shiny is installed and go to "app.R" directory
jid=$(SINGULARITYENV_port_num=$port_num \
	sbatch --time=$time \
	--output=run_overlap.out \
	--wrap "singularity exec \
	--bind $proj_dir:/mnt \
	--bind $img_dir/scripts:/scripts -B ~/.Xauthority \
	$img_dir/$img_name /bin/sh /scripts/app_simg.sbatch"| cut -f 4 -d' ')

echo -e "\n\nYou need to have x11 display server such as Xming running.\n"

node=$(squeue -j $jid -o "%R" -h)
status=$(squeue -j $jid -o "%t" -h)
while [[ $status != "R"* ]];do
	echo -e "Please wait setting up shiny app ....\n"
	sleep 10
	node=$(squeue -j $jid -o "%R" -h)
	status=$(squeue -j $jid -o "%t" -h)
done

echo Please ignore \"Failed to open connection..\" message.
echo ""
echo If you received \"connection to $node closed by remote host\",
echo please re-run run_overlap.sh.
echo ""
echo -e "After typing in your password, please wait until firefox appears ....\n"
echo -e "It may take some time for firefox to be responsive as it tries to load all the data\n"
echo -e "If a pop-up window appears, click on \"Create New Profile\"\n"
sleep 5
ssh -tX "$node" 'export port_num='"'$port_num'"'; \
        export img_dir='"'$img_dir'"'; \
        export img_name='"'$img_name'"'; \
        bash $img_dir/scripts/run_firefox.sh'
scancel $jid



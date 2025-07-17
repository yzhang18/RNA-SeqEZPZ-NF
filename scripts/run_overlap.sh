#!/bin/bash

# script to run shiny app
# How to run
# cd <project_dir>
# bash /export/apps/opt/rnaseq-pipeline/2.0/scripts/run_overlap.sh
# Examples:
# bash /export/apps/opt/rnaseq-pipeline/2.0/scripts/run_overlap.sh
# <project_dir> is where the fastq and outputs directory are
# project directory where the fastq and outputs directory are
#
#

# clear python path to prevent mixed up of python packages
unset PYTHONPATH

proj_dir=$(pwd)/../..
work_dir=$(pwd)

# singularity image directory
# find based on location of this script
img_dir=$HOME/singularity

# singularity image name
img_name=olechwin-pipelines-rnaseq.img

# IMPORTANT: It is assumed that:
# scripts to run analysis are in $img_dir/scripts
# reference to run analysis are in $img_dir/ref

echo -e "\nUsing singularity image and scripts in:" ${img_dir} "\n"
echo ""

# find available port
FROM=1035
TO=65535
port_num=$(singularity exec $img_dir/$img_name comm -23 \
<(seq "$FROM" "$TO" | sort) \
<(ss -Htan | awk '{print $4}' | cut -d':' -f2 | sort -u) \
| shuf | head -n 1)


APPTAINERENV_port_num=$port_num singularity exec \
    --bind $proj_dir:/mnt \
    --bind $work_dir:/scripts -B ~/.Xauthority \
    $img_dir/$img_name /bin/sh /scripts/app_simg.sh &

pid=$(ps aux | grep [s]hiny | awk '{print $2}')

# Activate environment where shiny is installed and go to "app.R" directory


echo Please ignore \"Failed to open connection..\" message.
echo ""
echo If you received \"connection to $node closed by remote host\",
echo please re-run run_overlap.sh.
echo ""
echo -e "After typing in your password, please wait until firefox appears ....\n"
echo -e "It may take some time for firefox to be responsive as it tries to load all the data\n"
echo -e "If a pop-up window appears, click on \"Create New Profile\"\n"
sleep 5

APPTAINERENV_port_num=$port_num singularity exec \
    $img_dir/$img_name bash -c "source activate firefox_env; \
    firefox --no-remote --new-window -P \'default\' http://127.0.0.1:$port_num"

kill -n 9 $pid


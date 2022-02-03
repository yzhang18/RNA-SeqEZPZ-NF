#!/bin/bash

# script to run shiny app
# How to run
# cd <project_dir>
# bash /apps/opt/rnaseq-pipeline/scripts/run_overlap.sh
# Examples:
# bash /apps/opt/rnaseq-pipeline/scripts/run_overlap.sh
# <project_dir> is where the fastq and outputs directory are
#set -x
# project directory where the fastq and outputs directory are
proj_dir=$(pwd)
cd $proj_dir

work_dir=$proj_dir/outputs

# clear python path to avoid reading in user's package sites
unset PYTHONPATH

# singularity image directory
# find based on location of this script
img_dir=$(dirname $(dirname $(readlink -f $0)))

# singularity image name
img_name=rnaseq-pipe-container.sif
# IMPORTANT: It is assumed that:
# scripts to run analysis are in $img_dir/scripts
# reference to run analysis are in $img_dir/ref

# find available port
FROM=1035
TO=65535
port_num=$(singularity exec $img_dir/$img_name comm -23 \
<(seq "$FROM" "$TO" | sort) \
<(ss -Htan | awk '{print $4}' | cut -d':' -f2 | sort -u) \
| shuf | head -n 1)

# Activate environment where shiny is installed and go to "app.R" directory
jid=$(SINGULARITYENV_port_num=$port_num \
	sbatch --time=5:00:00 \
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
echo -e "If a pop-up window appears, click on \"Create New Profile\"\n"
sleep 5
ssh -tX "$node" 'export port_num='"'$port_num'"'; \
        export img_dir='"'$img_dir'"'; \
        export img_name='"'$img_name'"'; \
        bash $img_dir/scripts/run_firefox.sh'
scancel $jid

#!/bin/bash

# How to run
# cd <project_dir>
# bash /gpfs0/home1/gdlessnicklab/cxt050/Steve/virtual_server/rnaseq-singularity/scripts/run_overlap.sh &> run_overlap.out &
# Examples:
# bash /gpfs0/home1/gdlessnicklab/cxt050/Steve/virtual_server/rnaseq-singularity/scripts/run_overlap.sh &> run_overlap.out &
# <project_dir> is where the fastq and outputs directory are
#set -x
# project directory where the fastq and outputs directory are
proj_dir=$(pwd)
cd $proj_dir

work_dir=$proj_dir/outputs

# singularity image directory
img_dir=/gpfs0/home1/gdlessnicklab/cxt050/Steve/virtual_server/rnaseq-singularity

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
#re='r1pl-hpc'
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
ssh -Y "$node" SINGULARITYENV_port_num=$port_num \
	singularity exec $img_dir/$img_name bash -c "source activate firefox_env ;\
	firefox --no-remote --new-window -P \"default\" http://127.0.0.1:$port_num"
scancel $jid

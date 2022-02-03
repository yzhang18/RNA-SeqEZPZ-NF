#!/bin/bash

# script to create genome STAR index
# How to run:
# cd <project_dir>
# bash /apps/opt/rnaseq-pipeline/scripts/run_star_index.sh &> run_star_index.out &
# Examples:
# cd ~/project1
# bash /apps/opt/rnaseq-pipeline/scripts/run_star_index.sh &> run_star_index.out &
#
# or to do nothing but echo all commands:
# bash /apps/opt/rnaseq-pipeline/scripts/run_star_index.sh echo &> run_star_index.out &
#

# set -x
date
# set 'run' to echo to simply echoing all commands
# set to empty to run all commands
if [ "$#" -eq 1 ] ; then
  run=$1
else
  run=
fi
# project directory
proj_dir=$(pwd)
cd $proj_dir

# specify number of cpus for star, bamCompare and bigwigCompare
ncpus=20

# singularity image directory
# find based on location of this script
img_dir=$(dirname $(dirname $(readlink -f $0)))

# singularity image name
img_name=rnaseq-pipe-container.sif
# IMPORTANT: It is assumed that:
# scripts to run analysis are in $img_dir/scripts
# reference to run analysis are in $img_dir/ref

# specify version of reference genome options hg19 or hg38
ref_ver=hg38
if [[ $ref_ver == 'hg19' ]]; then
	gtf_file=Homo_sapiens.GRCh37.75.gtf
	fasta_file=human_g1k_v37_decoy.fasta
elif [[ $ref_ver == 'hg38' ]]; then
	fasta_file=hg38.analysisSet.fa
	gtf_file=hg38.refGene.gtf
fi
# this is where star index will be stored
work_dir=$img_dir/ref/${ref_ver}/STAR_index
echo "all outputs will be stored in $work_dir"
if [ ! -d $work_dir ]; then
	mkdir -p $work_dir
fi
log_dir=$work_dir/logs
if [ ! -d $log_dir ]; then
	mkdir -p $log_dir
fi
echo "See $proj_dir/run_star_index.out to check the analysis progress"
echo "all other logs will be stored in $log_dir"
echo "log files contain all the commands run"
echo ""

#### generate index ####

tmp_jid=$(SINGULARITYENV_PYTHONPATH= \
	SINGULARITYENV_run=$run \
		SINGULARITYENV_ncpus=$ncpus \
		SINGULARITYENV_prefix=$prefix \
		SINGULARITYENV_ref_ver=$ref_ver \
		SINGULARITYENV_fasta_file=$fasta_file \
		SINGULARITYENV_gtf_file=$gtf_file \
		$run sbatch --output=$log_dir/star_index.out \
			--cpus-per-task $ncpus \
			--partition=himem \
			--mail-type=FAIL \
			--mail-user=$email \
			--job-name=star \
			--wrap "singularity exec \
				--bind $proj_dir:/mnt \
				--bind $img_dir/scripts:/scripts \
				--bind $img_dir/ref:/ref \
				$img_dir/$img_name \
				/bin/bash /scripts/star_index_simg.sbatch"| cut -f 4 -d' ')
echo "Generating STAR index job id: $tmp_jid"
echo ""
#set +x	


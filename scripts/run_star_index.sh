#!/bin/bash

# How to run:
# cd <project_dir>
# bash /gpfs0/home1/gdlessnicklab/cxt050/Steve/virtual_server/rnaseq-singularity/scripts/run_star_index.sh &> run_star_index.out &
# Examples:
# cd ~/Steve/virtual_server/project1
# bash /gpfs0/home1/gdlessnicklab/cxt050/Steve/virtual_server/rnaseq-singularity/scripts/run_star_index.sh &> run_star_index.out &
# or to do nothing but echo all commands:
# bash /gpfs0/home1/gdlessnicklab/cxt050/Steve/virtual_server/rnaseq-singularity/scripts/run_star_index.sh echo &> run_star_index.out &

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
img_dir=/gpfs0/home1/gdlessnicklab/cxt050/Steve/virtual_server/rnaseq-singularity
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

tmp_jid=$(SINGULARITYENV_run=$run \
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


# n_rep=$(printf "%s\n" "${repname_array[@]}" | sort -u | wc -l)
# echo number of replicates $n_rep
# if [[ n_rep -gt 1 ]]; then
	# #echo There is replicates	
	# #set -x
	# cd $work_dir
	# # dummy sbatch waiting for previous jobs to finish
	# tmp=$($run sbatch --wait \
	# --dependency=afterok:$jid3 \
	# --output=$log_dir/dummy.txt \
	# --wrap "echo dummy sbatch after bamtobw is finished")
	# # initialize job ids
	# jid3b=
	# for file in *_${repname_array[0]}_spike.bw; do
		# # set -x
		# # combined the replicates of bigwig into one bigwig
		# # echo $file
		# rep_files_array=($(ls ${file/_${repname_array[0]}_spike.bw/*_spike.bw} |\
			# grep -v $spikename | grep -v comb_spike.bw ))
		# j=1;
		# all_files_rep_cmd=
		# for i in "${rep_files_array[@]}"; do
		   # all_files_rep_cmd=$(echo $all_files_rep_cmd )
		   # b=$(echo -b$j | tr -d '[:space:]')
		   # all_files_rep_cmd=$(echo $all_files_rep_cmd $b )
		   # all_files_rep_cmd=$(echo $all_files_rep_cmd $i )
		   # j=$((j+1))
		# done
		# all_files_rep_cmd="$all_files_rep_cmd"
		# out_prefix=${file/_${repname_array[0]}_spike.bw/}
		# cd $proj_dir
		# #	
		# tmp_jid=$(SINGULARITYENV_run=$run \
				# SINGULARITYENV_all_files_rep_cmd=$all_files_rep_cmd \
				# SINGULARITYENV_ncpus=$ncpus \
				# SINGULARITYENV_out_prefix=$out_prefix \
				# $run sbatch --output=$log_dir/combinebw_${out_prefix}.out \
					# --dependency=afterok:$jid3 \
					# --cpus-per-task $ncpus \
					# --partition=himem \
					# --mail-type=FAIL \
					# --mail-user=$email \
					# --job-name=combinebw \
					# --wrap "singularity exec \
						# --bind $proj_dir:/mnt \
						# --bind $img_dir/scripts:/scripts \
						# $img_dir/cut-n-tag.sif \
							# /bin/bash /scripts/combinebw_simg.sbatch"| cut -f 4 -d' ')
		# echo "Combining $out_prefix bigwig replicates job id: $tmp_jid"
		# echo ""
		# if [ -z "$jid3b" ]; then
			# jid3b=$tmp_jid
		# else
			# jid3b=${jid3b},${tmp_jid}
		# fi
		# cd $work_dir
		# set +x
	# done
	# tmp=$($run sbatch --wait \
		# --dependency=afterok:$jid3b \
		# --output=$log_dir/dummy.txt \
		# --mail-type=END \
		# --mail-user=$email \
		# --job-name=run_create_tracks \
		# --wrap "echo dummy sbatch after combinebw is finished")
# else
	# echo There are no replicates
	# tmp=($run sbatch --wait \
		# --dependency=afterok:$jid3 \
		# --output=$log_dir/dummy.txt \
		# --mail-type=END \
		# --mail-user=$email \
		# --job-name=run_create_tracks \
		# --wrap "echo dummy sbatch after bamtobw is finished. There are no replicates")
# fi
# echo "Done creating spike-in adjusted visualization tracks"
# echo ""
# echo Tracks are \"*spike.bw\" in $work_dir
# echo ""
# date

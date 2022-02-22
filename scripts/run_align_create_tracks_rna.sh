#!/bin/bash

# How to run:
# cd <project_dir>
# bash /export/apps/opt/rnaseq-pipeline/2.0/scripts/run_align_create_tracks_rna.sh &> run_align_create_tracks_rna.out &
# Examples:
# cd project1
# bash /export/apps/opt/rnaseq-pipeline/2.0/scripts/run_align_create_tracks_rna.sh &> run_align_create_tracks_rna.out &
#
# by default, alignment is done to human reference genome hg19 unless specified genome=hg38:
# bash /export/apps/opt/rnaseq-pipeline/2.0/scripts/run_align_create_tracks_rna.sh genome=hg38 &> run_align_create_tracks_rna.out &
#
# or to run with specific time limit:
# bash /export/apps/opt/rnaseq-pipeline/2.0/scripts/run_align_create_tracks_rna.sh time=DD-HH:MM:SS &> run_align_create_tracks_rna.out &
#
# or to do nothing but echo all commands:
# bash /export/apps/opt/rnaseq-pipeline/2.0/scripts/run_align_create_tracks_rna.sh run=echo &> run_align_create_tracks_rna.out &
#
# or to run and printing all trace commands (i.e. set -x):
# bash /export/apps/opt/rnaseq-pipeline/2.0/scripts/run_align_create_tracks_rna.sh run=debug &> run_align_create_tracks_rna.out &

#set -x
set -e

# clear python path to avoid reading in user's package sites
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

	if [[ $1 == "help" ]] ;then
		echo ""
		echo 'usage: bash /export/apps/opt/rnaseq-pipeline/2.0/scripts/run_align_create_tracks_rna.sh [OPTION] &> run_align_create_tracks_rna.out &'
		echo ''
		echo DESCRIPTION
		echo -e '\trun alignment and tracks creation'
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
                echo -e "\tif using genome other than hg19 or hg38, need to put .fa or .fasta and gtf files in"
		echo -e "\tref/<genome-name> dir and set genome=<genome-name>."
		echo -e "time=1-00:00:00"
		echo -e "\tset SLURM time limit time=DD-HH:MM:SS, where ‘DD’ is days, ‘HH’ is hours, etc."
		echo -e "\tDefault is 1 day\n"
		exit
	fi
done

date
# converting samples.txt to unix format to remove any invisible extra characters
dos2unix -k samples.txt &> /dev/null

# set default parameters
# parameter to turn set -x on inside scripts
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
        run=
        debug=1
fi

echo -e "\nAligning reads to $ref_ver and create tracks for visualization\n"
echo -e "Options used to run:"
echo padj="$padj"
echo time="$time"
echo genome="$ref_ver"
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
echo "See $proj_dir/run_align_create_tracks_rna.out to check the analysis progress"
echo "All other logs including scripts ran will be stored in $log_dir"
echo ""

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

# copying this script for records
$(cp $img_dir/scripts/run_align_create_tracks_rna.sh $log_dir/run_align_create_tracks_rna.sh)

### specify reference genome
# if ref genome is not in $img_dir/ref, set genome_dir to  project dir
genome_dir=$img_dir/ref/$ref_ver
if [[ ! -d $genome_dir ]];then
        genome_dir=$proj_dir/ref/$ref_ver
fi
echo $genome_dir
gtf_file=$(find $genome_dir -name *.gtf | xargs basename)
fasta_file=$(find $genome_dir -name *.fasta -o -name *.fa | xargs basename)
chr_info=$(find $genome_dir -name *.chrom.sizes | xargs basename)
star_index_dir=$genome_dir/STAR_index

# calculate genome_size
genome_size=$(grep -v ">" $genome_dir/$fasta_file | grep -v "N" | wc | awk '{print $3-$1}')
chr_info=$(find $genome_dir -name *.chrom.sizes | xargs basename)

# getting samples info from samples.txt
$(sed -e 's/[[:space:]]*$//' samples.txt | sed 's/"*$//g' | sed 's/^"*//g' > samples_tmp.txt)
$(mv samples_tmp.txt samples.txt)
groupname_array=($(awk '!/#/ {print $1}' samples.txt))
repname_array=($(awk '!/#/ {print $3}' samples.txt))
email=$(awk '!/#/ {print $5;exit}' samples.txt | tr -d '[:space:]')
filename_string_array=($(awk '!/#/ {print $6}' samples.txt))
string_pair1_array=($(awk '!/#/ {print $7}' samples.txt))
string_pair2_array=($(awk '!/#/ {print $8}' samples.txt))

#### alignment ####
# STAR first pass
# initialize job ids
jid3=
for i in "${!groupname_array[@]}"; do 
	groupname=${groupname_array[$i]}
	repname=${repname_array[$i]}
	string_pair1=${string_pair1_array[$i]}
	string_pair2=${string_pair2_array[$i]}
	filename_string=${filename_string_array[$i]}
	prefix=${groupname}_${repname}
	# get the files for the same groupname and replicate
	# (i.e from multiple lanes)
	cd $work_dir/trim
	file=($(ls *val_1.fq.gz | \
		awk -v groupname=$groupname -v repname=$repname -v filename_string=$filename_string \
		'$1 ~ groupname && $1 ~ repname && $1 ~ filename_string {print $1}'))
	# combining files from different lanes
	read1=$(printf ",%s" "${file[@]}")
	read1=${read1:1}
	file=($(ls *val_2.fq.gz | \
		awk -v groupname=$groupname -v repname=$repname -v filename_string=$filename_string \
		'$1 ~ groupname && $1 ~ repname && $1 ~ filename_string {print $1}'))
	read2=$(printf ",%s" "${file[@]}")
	read2=${read2:1}
	cd $proj_dir
	#echo $read1
	#echo $read2
	cd $proj_dir
	# set -x
	tmp_jid=$(SINGULARITYENV_PYTHONPATH= \
	SINGULARITYENV_read1=$read1 \
		SINGULARITYENV_read2=$read2 \
		SINGULARITYENV_run=$run \
		SINGULARITYENV_ncpus=$ncpus \
		SINGULARITYENV_prefix=$prefix \
		SINGULARITYENV_ref_ver=$ref_ver \
		$run sbatch --output=$log_dir/star_pass1_${prefix}.out \
			--cpus-per-task $ncpus \
			--partition=himem \
			--mail-type=FAIL \
			--mail-user=$email \
			--job-name=star_pass1 \
			--time=$time \
			--wrap "singularity exec \
				--bind $proj_dir:/mnt \
				--bind $img_dir/scripts:/scripts \
				--bind $genome_dir:/ref \
				$img_dir/$img_name \
				/bin/bash /scripts/star_pass1_simg.sbatch"| cut -f 4 -d' ')
	echo "Running STAR first pass for $prefix job id: $tmp_jid"
	echo ""
	if [ -z "$jid3" ]; then
		jid3=$tmp_jid
	else
		jid3=${jid3}:${tmp_jid}
	fi
#set +x
done

# dummy job to wait for STAR first pass to finish before
# continuing on
tmp=$($run sbatch \
	--dependency=afterok:$jid3 \
	--output=$log_dir/dummy.txt \
	--job-name=after_star_pass1 \
	--time=5:00 \
	--wrap "echo dummy sbatch to wait for STAR 1st pass to be finished"| cut -f 4 -d' ')

# message if jobs never satisfied or cancelled
check_jid3=$(echo $jid3 | sed 's/:/,/g')
state=($(squeue -j $check_jid3 -h))

while [ ${#state[@]} -ne 0 ];
do
        sleep 10
        state=($(squeue -j $check_jid3 -h))
done

reason=$(squeue -j $tmp -o "%R" -h)
state=$(sacct -j $tmp --format=state | tail -n +3 | head -n 1)
if [[ $reason == *"DependencyNeverSatisfied"* || $state == *"CANCELLED"* ]]; then
	scancel $tmp
	echo -e "STAR 1st pass run failed. Please check star_pass1_* files in $log_dir\n"
	exit
fi
# STAR second pass
genomeForPass2=$work_dir/STAR_2pass/GenomeForPass2

# get all the junctions files for input to 2nd pass alignment
cd $genomeForPass2
sjFiles=($(find `pwd` -name '*_SJ.out.tab.Pass1.sjdb'))
sjFiles=$(printf " %s" "${sjFiles[@]}")
sjFiles=${sjFiles:1}
sjFiles=${sjFiles//$genomeForPass2/..\/STAR_2pass\/GenomeForPass2}

cd $proj_dir
# initialize job ids
jid4=
for i in "${!groupname_array[@]}"; do
	groupname=${groupname_array[$i]}
	repname=${repname_array[$i]}
	string_pair1=${string_pair1_array[$i]}
	string_pair2=${string_pair2_array[$i]}
	filename_string=${filename_string_array[$i]}
	prefix=${groupname}_${repname}
	# get the files for the same groupname and replicate
	# (i.e from multiple lanes)
	cd $work_dir/trim
	file=($(ls *val_1.fq.gz | \
		awk -v groupname=$groupname -v repname=$repname -v filename_string=$filename_string \
		'$1 ~ groupname && $1 ~ repname && $1 ~ filename_string {print $1}'))
	read1=$(printf ",%s" "${file[@]}")
	read1=${read1:1}
	file=($(ls *val_2.fq.gz | \
		awk -v groupname=$groupname -v repname=$repname -v filename_string=$filename_string \
		'$1 ~ groupname && $1 ~ repname && $1 ~ filename_string {print $1}'))
	read2=$(printf ",%s" "${file[@]}")
	read2=${read2:1}
	cd $proj_dir
	#echo $read1
	#echo $read2
	# set -x
	tmp_jid=$(SINGULARITYENV_PYTHONPATH= \
	SINGULARITYENV_read1=$read1 \
			SINGULARITYENV_read2=$read2 \
			SINGULARITYENV_run=$run \
			SINGULARITYENV_ncpus=$ncpus \
			SINGULARITYENV_sjFiles=$sjFiles \
			SINGULARITYENV_prefix=$prefix \
			SINGULARITYENV_ref_ver=$ref_ver \
			SINGULARITYENV_fasta_file=$fasta_file \
			SINGULARITYENV_gtf_file=$gtf_file \
			SINGULARITYENV_genome_size=$genome_size \
			$run sbatch --output=$log_dir/star_pass2_${prefix}.out \
				--cpus-per-task $ncpus \
				--partition=himem \
				--mail-type=FAIL \
				--dependency=afterok:$jid3 \
				--mail-user=$email \
				--job-name=star_pass2 \
				--mem=128G \
				--time=$time \
				--wrap "singularity exec \
					--bind $proj_dir:/mnt \
					--bind $img_dir/scripts:/scripts \
					--bind $genome_dir:/ref \
					$img_dir/$img_name \
					/bin/bash /scripts/star_pass2_simg.sbatch"| cut -f 4 -d' ')
	echo "Running STAR second pass for $prefix job id: $tmp_jid"
	echo ""
	if [ -z "$jid4" ]; then
		jid4=$tmp_jid
	else
		jid4=${jid4}:${tmp_jid}
	fi
	#set +x	
done	

# if there are replicates, combine the bw files for each replicate 
n_rep=$(printf "%s\n" "${repname_array[@]}" | sort -u | wc -l)
echo number of replicates $n_rep
if [[ n_rep -gt 1 ]]; then
	#echo There is replicates	
	#set -x
	cd $work_dir
	# dummy sbatch waiting for previous jobs to finish
	tmp=$($run sbatch \
	--dependency=afterok:$jid4 \
	--output=$log_dir/dummy.txt \
	--time=5:00 \
	--wrap "echo dummy sbatch waiting for STAR 2nd pass to be finished"| cut -f 4 -d' ')
	
	# message if jobs never satisfied or cancelled
	check_jid4=$(echo $jid4 | sed 's/:/,/g')
	state=($(squeue -j $check_jid4 -h))

	while [ ${#state[@]} -ne 0 ];
	do
        	sleep 10
        	state=($(squeue -j $check_jid4 -h))
	done
				
	reason=$(squeue -j $tmp -o "%R" -h)
	state=$(sacct -j $tmp --format=state | tail -n +3 | head -n 1)
	if [[ $reason == *"DependencyNeverSatisfied"* || $state == *"CANCELLED"* ]]; then
		scancel $tmp
		echo -e "STAR 2nd pass run failed. Please check star_pass2_* files in $log_dir\n"
		exit
	fi	
	
	# initialize job ids
	jid4b=
	cd $work_dir/bw_files
	for file in *_${repname_array[0]}_norm.bw; do
		# set -x
		# combined the replicates of bigwig into one bigwig
		# echo $file
		rep_files_array=($(ls ${file/_${repname_array[0]}_norm.bw/*_norm.bw}))
		all_files_rep_cmd=$(printf " %s" "${rep_files_array[@]}")
		all_files_rep_cmd=${all_files_rep_cmd:1}
		echo $all_files_rep_cmd
		# find groupname
		cd $proj_dir
		#echo $chr_info
		groupname=$(awk -v file=$file 'file ~ $1 {print $1;exit}' samples.txt)
		tmp_jid=$(SINGULARITYENV_PYTHONPATH= \
		SINGULARITYENV_run=$run \
				SINGULARITYENV_all_files_rep_cmd=$all_files_rep_cmd \
				SINGULARITYENV_ncpus=$ncpus \
				SINGULARITYENV_groupname=$groupname \
				SINGULARITYENV_ref_ver=$ref_ver \
				SINGULARITYENV_fasta_file=$fasta_file \
				SINGULARITYENV_chr_info=$chr_info \
				$run sbatch --output=$log_dir/combinebw_${groupname}.out \
					--dependency=afterok:$jid4 \
					--cpus-per-task $ncpus \
					--partition=himem \
					--mail-type=FAIL \
					--mail-user=$email \
					--job-name=combinebw \
					--time=$time \
					--wrap "singularity exec \
						--bind $proj_dir:/mnt \
						--bind $genome_dir:/ref \
						--bind $img_dir/scripts:/scripts \
						$img_dir/$img_name \
							/bin/bash /scripts/combinebw_simg.sbatch"| cut -f 4 -d' ')
		echo "Combining $out_prefix bigwig replicates job id: $tmp_jid"
		echo ""
		cd $work_dir/bw_files
		if [ -z "$jid4b" ]; then
			jid4b=$tmp_jid
		else
			jid4b=${jid4b}:${tmp_jid}
		fi
		# set +x
	done
	jid4c=$($run sbatch \
		--dependency=afterok:$jid4b \
		--output=$log_dir/dummy.txt \
		--mail-type=END \
		--mail-user=$email \
		--time=5:00 \
		--job-name=run_align_create_tracks_rna \
		--wrap "echo dummy sbatch after combinebw is finished"| cut -f 4 -d' ')
	# message if jobs never satisfied or canceled
	check_jid4b=$(echo $jid4b | sed 's/:/,/g')
	state=($(squeue -j $check_jid4b -h))

	while [ ${#state[@]} -ne 0 ];
	do
        	sleep 10
        	state=($(squeue -j $check_jid4b -h))
	done
				
	reason=$(squeue -j $jid4c -o "%R" -h)
	state=$(sacct -j $jid4c --format=state | tail -n +3 | head -n 1)
	if [[ $reason == *"DependencyNeverSatisfied"* || $state == *"CANCELLED"* ]]; then
		scancel $jid4c
		echo -e "Combining bw files failed. Please check combinebw_* files in $log_dir\n"
		exit
	fi	
else
	echo There are no replicates
	jid4c=$($run sbatch \
		--dependency=afterok:$jid4 \
		--output=$log_dir/dummy.txt \
		--mail-type=END \
		--mail-user=$email \
		--time=5:00 \
		--job-name=run_align_create_tracks_rna \
		--wrap "echo dummy sbatch after STAR 2nd pass is finished. There are no replicates"| cut -f 4 -d' ')
	# message if jobs never satisfied or cancelled
	check_jid4=$(echo $jid4 | sed 's/:/,/g')
	state=($(squeue -j $check_jid4 -h))

	while [ ${#state[@]} -ne 0 ];
	do
        	sleep 10
        	state=($(squeue -j $check_jid4 -h))
	done
			
	reason=$(squeue -j $jid4c -o "%R" -h)
	state=$(sacct -j $jid4c --format=state | tail -n +3 | head -n 1)
	if [[ $reason == *"DependencyNeverSatisfied"* || $state == *"CANCELLED"* ]]; then
		scancel $jid4c
		echo -e "STAR 2nd pass run failed. Please check star_pass2_* files in $log_dir\n"
		exit
	fi	
fi

message="Done alignment using STAR 2-pass approach and created bw files for visualization\n\n\
Tracks for each replicate and combined replicates are \"*norm.bw\" and \"comb.bw\" respectively\n\
in $work_dir\n\
Aligned bam files are in $work_dir/STAR_2pass/Pass2\n\n"

tmp=$($run sbatch --dependency=afterok:$jid4c \
		--output=$log_dir/dummy.txt \
		--mail-type=END \
		--mail-user=$email \
		--time=5:00 \
		--job-name=run_align_create_tracks_rna \
		--export message="$message",proj_dir=$proj_dir \
		--wrap "echo -e \"$message\"$(date) >> $proj_dir/run_align_create_tracks_rna.out"| cut -f 4 -d' ')


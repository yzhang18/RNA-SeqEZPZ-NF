#!/bin/bash

set -e
#set -v
# This is for running the second pass of STAR alignment
echo Hello from $SINGULARITY_CONTAINER
cd /mnt/outputs/trim
source activate rna_env
echo "STAR version " $(STAR --version)
set -x
genomeDir=/ref/STAR_index
mkdir -p ../STAR_2pass/Pass2
runDir=../STAR_2pass/Pass2
export OMP_NUM_THREADS=$ncpus
echo $read1 $read2
$run STAR --genomeDir $genomeDir \
	--readFilesIn $read1 $read2 \
	--runThreadN $ncpus \
	--readFilesCommand zcat \
	--outSAMtype BAM SortedByCoordinate Unsorted \
	--sjdbFileChrStartEnd $sjFiles \
	--limitSjdbInsertNsj 4000000 \
	--outFileNamePrefix "$runDir/${prefix}_"
set +x
conda deactivate

# create bam index and normalized bams
source activate samtools_env
samtools --version
cd $runDir
runDir=$(pwd)
outdir=/mnt/outputs/bw_files
mkdir -p $outdir
set -x
$run samtools index ${prefix}_Aligned.sortedByCoord.out.bam
set +x
conda deactivate

# calculating effective genome size
source activate khmer_env
set -x
cd /opt/khmer
unique-kmers.py -k $readlength $fasta_file &> /mnt/outputs/kmers.txt
genome_size=$(tail -n1 /mnt/outputs/kmers.txt | cut -f2 -d:)
set +x
conda deactivate

cd $runDir
source activate deeptools_env
bamCoverage --version
set -x
$run bamCoverage -b ${prefix}_Aligned.sortedByCoord.out.bam \
	-o $outdir/${prefix}_norm.bw \
	-p $ncpus \
	--binSize 10 \
	--smoothLength 20 \
	--normalizeUsing BPM \
	--effectiveGenomeSize $genome_size
set +x
conda deactivate

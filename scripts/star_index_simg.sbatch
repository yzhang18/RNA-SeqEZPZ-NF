#!/bin/bash

set -e
#set -v
echo Hello from $SINGULARITY_CONTAINER
source activate rna_env
echo "STAR version " $(STAR --version)
set -x
export OMP_NUM_THREADS=$ncpus
$run STAR --runMode genomeGenerate \
	--genomeDir /ref/STAR_index \
	--runThreadN $ncpus \
	--sjdbGTFfile /ref/$gtf_file \
	--genomeFastaFiles /ref/$fasta_file \
	--sjdbOverhang 100
set +x
conda deactivate
# create chromosome info file
source activate samtools_env
samtools --version
set -x
chr_info="${fasta_file%.*}"
$run samtools faidx /ref/$fasta_file
cut -f1,2 /ref/${fasta_file}.fai > /ref/${chr_info}.chrom.sizes
set +x
conda deactivate

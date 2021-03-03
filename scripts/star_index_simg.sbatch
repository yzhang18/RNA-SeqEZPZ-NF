#!/bin/bash

set -e
#set -v
echo Hello from $SINGULARITY_CONTAINER
source activate rna_env
echo "STAR version " $(STAR --version)
set -x
genomeDir=/ref/$ref_ver
export OMP_NUM_THREADS=$ncpus
STAR --runMode genomeGenerate \
	--genomeDir $genomeDir/STAR_index \
	--runThreadN $ncpus \
	--sjdbGTFfile $genomeDir/$gtf_file \
	--genomeFastaFiles $genomeDir/$fasta_file \
	--sjdbOverhang 100 
set +x
conda deactivate
	
	

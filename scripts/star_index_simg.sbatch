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
	--sjdbGTFfile /ref_gtf/$target_gtf_name \
	--genomeFastaFiles /ref_fa/$target_fa_name \
	--sjdbOverhang 100
set +x
conda deactivate
# create chromosome info file
source activate samtools_env
samtools --version
set -x
chr_info="${fasta_file%.*}"
$run samtools faidx /ref_fa/$target_fa_name
cut -f1,2 /ref_fa/${target_fa_name}.fai > /ref/${chr_info}.chrom.sizes
set +x
conda deactivate

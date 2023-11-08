#!/bin/bash

set -e
#set -v
echo Hello from $SINGULARITY_CONTAINER
source activate rna_env
echo "STAR version " $(STAR --version)
set -x
export OMP_NUM_THREADS=$ncpus
cd /mnt
# only generate STAR index if it doesn't exist yet
if [[ ! "$(ls -A /ref/STAR_index)" ]]; then
	$run STAR --runMode genomeGenerate \
		--genomeDir /ref/STAR_index \
		--runThreadN $ncpus \
		--sjdbGTFfile /ref_gtf/$target_gtf_name \
		--genomeFastaFiles /ref_fa/$target_fa_name \
		--sjdbOverhang 100
fi
set +x
conda deactivate
# create chromosome info file
source activate samtools_env
samtools --version
set -x
chr_info="${target_fa_name%.*}"
# check if fai file exist
if [[ ! -f /ref_fa/${target_fa_name}.fai ]]; then
	$run samtools faidx /ref_fa/$target_fa_name
fi
# check if chrom.sizes exist
if [[ ! -f /ref/${chr_info}.chrom.sizes ]];then
	cut -f1,2 /ref_fa/${target_fa_name}.fai > /ref/${chr_info}.chrom.sizes
fi
set +x
conda deactivate

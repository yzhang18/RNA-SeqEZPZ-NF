#!/bin/bash

#set -x
set -e
#set -v
echo Hello from $SINGULARITY_CONTAINER
echo singularity image version: $version

source activate rna_env
trim_galore --version
set -x
if [ ! -d /mnt/outputs/fastqc_rslt ]; then
	mkdir -p /mnt/outputs/fastqc_rslt
fi
if [ ! -d /mnt/outputs/trim ]; then
	mkdir -p /mnt/outputs/trim
fi
cd /mnt/outputs/trim
$run trim_galore -q 20 --length 20 \
	--illumina \
	--cores $ncpus_trim \
	--paired \
	$path_tech_r1 $path_tech_r2 \
	--basename $prefix \
	--fastqc --fastqc_args "--outdir /mnt/outputs/fastqc_rslt"
# 4 is thought to be optimal number of cores in trim_galore
set +x
conda deactivate

#!/bin/bash

. /gpfs0/home/gdlessnicklab/cxt050/opt/miniconda3/etc/profile.d/conda.sh
conda activate local_denovo_asm_hifi_env
set -x
cd ~/../lab/data/210913_Lessnick_GSL-JC-2349_fusions_mousecells
out_dir=~/Steve/virtual_server/rnaseq-singularity/project_ex3
mkdir -p $out_dir/fastq
for file in fastq/*.gz;do seqtk sample -s100 $file 10000|gzip > ${out_dir}/${file/fastq.gz/sub.fastq.gz} ;done
set +x
conda deactivate

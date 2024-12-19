#!/bin/bash

set -e
#set -v
# set -x
echo Hello from $SINGULARITY_CONTAINER
cd /mnt/outputs
if [ ! -d fastqc_rslt ]; then
	mkdir fastqc_rslt
fi
cd fastqc_rslt
source activate rna_env
multiqc --version
set -x
$run multiqc -f $input_dir
set +x
conda deactivate

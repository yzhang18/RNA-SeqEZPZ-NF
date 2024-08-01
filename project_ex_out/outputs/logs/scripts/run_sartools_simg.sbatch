#!/bin/bash

# script to run template_script_DESeq2_simg.R

set -e
echo hello from $SINGULARITY_CONTAINER
source activate sartools_env
set -x
cd /mnt/outputs
$run Rscript /scripts/template_script_DESeq2_simg.R padj=$padj email=$email batch_adjust=$batch_adjust
set +x
conda deactivate
#!/bin/bash

set -e
echo Hello from $SINGULARITY_CONTAINER
cd /mnt/outputs/bw_files
source activate wiggletools_env
set -x
$run wiggletools mean $all_files_rep_cmd | \
	wigToBigWig stdin /ref/$ref_ver/$chr_info ${groupname}_comb.bw
set +x
conda deactivate
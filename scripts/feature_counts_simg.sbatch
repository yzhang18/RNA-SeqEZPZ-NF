#!/bin/bash

set -e

echo Hello from $SINGULARITY_CONTAINER
set -x
cd /mnt/outputs/STAR_2pass/Pass2
set +x
source activate rna_env
featureCounts -v
set -x
$run featureCounts -p \
	-a $gtf_file \
	-t exon \
	-g gene_id \
	-T $ncpus \
	-o ${prefix}_counts.txt \
	${prefix}_Aligned.sortedByCoord.out.bam
set +x
conda deactivate

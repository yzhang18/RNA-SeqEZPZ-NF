#!/bin/bash

set -e
set -v

echo Hello from $SINGULARITY_CONTAINER
cd /mnt/outputs/STAR_2pass/Pass2
source activate rna_env
featureCounts -v
set -x
$run featureCounts -p \
	-a /ref_gtf/$target_gtf_name \
	-t exon \
	-g gene_id \
	-T $ncpus \
	-o ${prefix}_counts.txt \
	${prefix}_Aligned.sortedByCoord.out.bam
set +x
conda deactivate

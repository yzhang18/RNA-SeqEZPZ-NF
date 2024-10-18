#!/bin/bash

set -e

echo Hello from $SINGULARITY_CONTAINER
set -x
cd /mnt/outputs/STAR_2pass/Pass2
set +x
source activate rna_env
featureCounts -v
set -x
# --countReadPairs read-pairs will be counted instead of reads
# --countReadPairs needed for featureCounts v 2.0.6
# --countReadPairs NOT needed when using featureCounts v 2.0.1
$run featureCounts -p \
	--countReadPairs \
	-a $gtf_file \
	-t exon \
	-g gene_id \
	-T $ncpus \
	-o ${prefix}_counts.txt \
	${prefix}_Aligned.sortedByCoord.out.bam
set +x
conda deactivate

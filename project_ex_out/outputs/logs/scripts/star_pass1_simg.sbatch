#!/bin/bash

set -e
#set -v
echo Hello from $SINGULARITY_CONTAINER
cd /mnt/outputs/trim
source activate rna_env
echo "STAR version " $(STAR --version)
set -x
# binding /star_index_dir instead
#genomeDir=/ref/STAR_index
mkdir -p ../STAR_2pass/Pass1
runDir=../STAR_2pass/Pass1
export OMP_NUM_THREADS=$ncpus
$run STAR --genomeDir /star_index_dir \
	--readFilesIn $read1 $read2 \
	--runThreadN $ncpus \
	--readFilesCommand zcat \
	--outFileNamePrefix "$runDir/${prefix}_"

# preparing for STAR second pass
genomeForPass2=../STAR_2pass/GenomeForPass2
mkdir -p $genomeForPass2
# filtering non-canonical chromosomes
$run awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} {if($5>0){print $1,$2,$3,strChar[$4]}}' \
../STAR_2pass/Pass1/${prefix}_SJ.out.tab > $genomeForPass2/${prefix}_SJ.out.tab.Pass1.sjdb
set +x
conda deactivate

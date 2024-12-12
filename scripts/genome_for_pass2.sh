#!/bin/bash


pair_id=$1

awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} {if($5>0){print $1,$2,$3,strChar[$4]}}' ${pair_id}_SJ.out.tab > ${pair_id}_SJ.out.tab.Pass1.sjdb



#!/bin/bash
# $1: input tfbs
# $2: input count matrix
# $3: output file

awk 'NR>1' $2 > ${2}.tmp 
bedtools intersect -a $1 -b ${2}.tmp -wa -wb > $3

# clean up 
rm ${2}.tmp 


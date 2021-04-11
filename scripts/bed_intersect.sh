#!/bin/bash
# $1: aligned bed gz file
# $2: input tfbs file
# $3: output file
zcat $1 | awk -F '\t' 'BEGIN {OFS = FS} {$5=$1":"$2"-"$3;$2=int(($2+$3)/2); $3=$2; print}' |  bedtools intersect -a $2 -b - -wa -wb | gzip - > $3


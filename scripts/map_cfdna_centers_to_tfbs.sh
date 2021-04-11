#!/bin/bash
# $1 : bed file
# $2 : cfdna bed
# $3 : intersection and fragment length count file
# $4 : verbose file
# $5 : center map bedgz 

zcat $2 | awk '{OFS="\t"}{$2 = $2 - 1; print $0}' | bedtools intersect -a $1 -b - -sorted -wa -wb  | gzip - > ${5}

zcat ${5} | python scripts/create_cfdna_flen_list.py > $3 
zcat ${5} | cut -f1-6,13- | gzip - > $4 

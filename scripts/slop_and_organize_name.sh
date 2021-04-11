#!/bin/bash
# $1: input bed file
# $2: genome size file
# $3: slop value
# $4: output bed file

bedtools slop -b $3 -g $2 -i $1 > ${4}.tmp

awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"|"$4}' ${4}.tmp > $4  

# cleanup 
rm ${4}.tmp

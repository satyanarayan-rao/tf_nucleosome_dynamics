#!/bin/bash
# $1: bed file to order
# $2: output bed file 

awk '{if ($NF == 0) {print $0}}' $1 | tr '-' ' ' | sort -k5,5nr  > ${1}.tmp.0.bed
grep -w "inf" $1   > ${1}.tmp.inf.bed
awk '{if ($NF == 1) {print $0}}' $1  | grep -vw "inf" | tr '-' ' ' | sort -k5,5nr  > ${1}.tmp.1.bed

cat ${1}.tmp.inf.bed ${1}.tmp.1.bed ${1}.tmp.0.bed > ${2}.tmp.bed 
# get maximum distance of peak so that I can create a bed file in matrix form by slopping from center

max_dist=`awk '{print $5}' ${2}.tmp.bed | tr '-' ' ' | grep -v "inf" | sort -nr | head -1`
next_thousand=`echo "((${max_dist}/1000) + 1 )*1000" | bc` 

# create a mid_point_bed  
awk '{print $1"\t"int(($2 + $3)/2)"\t"int(($2 + $3)/2) +1"\t"$4"\t"$5"\t"$6}'  ${2}.tmp.bed > ${2}.mid_point.bed 
bedtools slop -b $next_thousand -g metadata/hg38.genome  -i ${2}.mid_point.bed > $2 

# cleanup 
#rm ${1}.tmp.inf.bed ${1}.tmp.1.bed ${1}.tmp.0.bed 
#rm ${2}.tmp.bed ${2}.mid_point.bed  


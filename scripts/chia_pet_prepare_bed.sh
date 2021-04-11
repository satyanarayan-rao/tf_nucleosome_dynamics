#!/bin/bash
# $1: head file 
# $2: tail file
# $3: output file
awk -F'[-:]' ' {OFS="\t"} {print $1, $2, $3}' $1 > ${1}.tmp 
paste ${1}.tmp ${2} | awk '{print $0"|head_tail"}' > ${1}.tmp1 

awk -F'[-:]' ' {OFS="\t"} {print $1, $2, $3}' $2 > ${2}.tmp 
paste ${2}.tmp ${1} | awk '{print $0"|tail_head"}' > ${2}.tmp1 

cat ${1}.tmp1 ${2}.tmp1 > ${3} 
rm ${1}.tmp ${1}.tmp1 ${2}.tmp ${2}.tmp1 

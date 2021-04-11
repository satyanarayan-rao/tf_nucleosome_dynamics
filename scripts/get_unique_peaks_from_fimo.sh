#!/bin/bash
# $1: fimo bed file
# $2: append lines to the count tsv file
cnt=`awk -F '[\t|]' '{print $4}' $1 | sort | uniq | wc -l`
cp ${2} ${2}.tmp
echo -e "$cnt ${1}_unique_peaks.bed" >> ${2}.tmp
sort -k1,1gr ${2}.tmp > ${2}

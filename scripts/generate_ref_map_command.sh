#!/bin/bash 
grep -v "^#" ${1} | sed '/^$/d' > ${1}.selected 
while read tfbs_bed fragment_bed len_col out_file chunk verbose
do
    echo Rscript scripts/cfdna_center_in_tfbs.R ${tfbs_bed} ${fragment_bed}  ${out_file} ${chunk} ${len_col} ${verbose}
done < ${1}.selected

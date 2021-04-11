#!/bin/bash
# $1: raw count matrix: tcga 
# $2: Patient id
# $3: output matrix
#

# first get the column ids of patients
fname_raw_count=`echo $1 | awk -F'/' '{print $NF}'`
head -1 $1 | tr '\t' '\n' > tmp/header_${fname_raw_count}
for i in `cat $2`; do grep -n "$i" tmp/header_${fname_raw_count} | \
    awk -F':' '{print $1}'; done > tmp/cl_id_${fname_raw_count}

fields_to_print=`paste -d',' -s tmp/cl_id_${fname_raw_count}`
first_five="1,2,3,4,5,"
all_selected_cols=`echo ${first_five}${fields_to_print}`
cut -f${all_selected_cols}  $1 > $3

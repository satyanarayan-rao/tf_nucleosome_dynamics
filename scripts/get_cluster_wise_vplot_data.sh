#!/bin/bash
# $1: cluster tagged intersection file
# $2: cl_id
# $3: output file
# $4: start bed file input_bed/50b_* : just to get column numbers
# $5: bedgz length column

zcat $1 | awk '{if ($1==cl_id) {print $0}}' cl_id=$2 | cut -f2- | gzip - > ${3}.tmp.gz

# get number of columns in bed
bed_cols=`awk '{print NF}' $4 | sort | uniq -c | awk '{print $2}'`
total=`echo "$bed_cols + $5" | bc -l`
python scripts/gen_vplot_data.py ${3}.tmp.gz $total $bed_cols $3 

rm ${3}.tmp.gz

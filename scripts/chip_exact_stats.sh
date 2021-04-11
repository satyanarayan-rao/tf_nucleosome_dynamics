#!/bin/bash
# input bed
# string to search in bed
# bw file 
# output exact stats
bname=`basename $1`
grep -w $2 $1 > tmp/${bname}.$2
python $NGS_SCRIPTS_DIR/map_bw_exact_stats_to_bed.py $3 tmp/${bname}.$2 $4

# clean up 
rm tmp/${bname}.$2 


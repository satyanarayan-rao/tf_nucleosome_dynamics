#!/bin/bash
# $1: chip bed
# $2: output shuf line file
# $3: output bed file
# $4: number of sites to be selected
num_lines=`wc -l $1 | awk '{print $1}'`
python scripts/shuffle_line_ids.py $num_lines $2
shuf --random-source=$2 -n $4 $1 > $3  

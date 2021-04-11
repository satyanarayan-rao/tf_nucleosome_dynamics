#!/bin/bash
# $1: input file (first four columns of bed and then min_len to max_len)
# $2: k: to ignore first k columns 
# $3: output file

start_col=`echo "$2 + 1" | bc`
echo $start_col

awk '{sum=0;for (i=start_col; i<=NF; i++){sum+=$i}; print $1"^"$2"^"$3"^"$4"\t"sum}' start_col=$start_col $1 > $3

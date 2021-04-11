#!/bin/bash
# $1: original vplot matrix for a cluster
# $2: output: zeros added on top till the first fargment length 

# first get the minimum fragment length

min_len=`head -1 $1 | awk '{print $1}'`

min_len_minus_one=`echo "$min_len - 1" | bc`
# get total fields  
xra=`head -1 $1 | awk '{print NF -1}'`

for i in `seq 0 ${min_len_minus_one}`
do 
    echo $i
    for j in `seq 1 $xra`
    do 
        echo 0
    done > 
  
done 

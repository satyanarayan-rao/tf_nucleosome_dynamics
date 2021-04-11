#!/bin/bash 

less_one=`echo "$# - 1" | bc `
last="$#"
echo $less_one
for i in  `seq $less_one`
do
    wc -l "${!i}" 

done > ${!last}

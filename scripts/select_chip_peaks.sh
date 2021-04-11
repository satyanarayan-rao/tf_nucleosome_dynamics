#!/bin/bash 
# $1: fimo bed - contains motif and peak information
# $2: boundary limits 
# $3: output bed file

awk '{print $4}' $1 | awk -F'[-:|]' '{print $1"\t"$2"\t"$3"\t"$4}' | sort -k1,1 -k2,2n -k3,3 --stable --uniq | awk '{if ( ($3 - $2) <=w)  {print $0}}' w=$2 > $3

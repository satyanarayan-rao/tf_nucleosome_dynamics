#!/bin/bash
# $1: starting bed file from public/inhouse sources
# $2: output bed file
# $3: bed annotation
# $4: part of the name field annotation

# head -5 /beevol/home/zukowski/ER_CUTnRUN/07302020_iterative_peaks/06.MCF7_E2_FOXA1_peaks.bed 
# 11	469876	469906	11_469891	625.099193906176
# 11	537456	537486	11_537471	625.099193906176
# 11	844316	844346	11_844331	625.099193906176
# 11	955606	955636	11_955621	357.199539374958
# 11	1220016	1220046	11_1220031	357.199539374958

check_chr=`head $1 | grep "^chr"`
if [ $? -eq 1 ]
then
    awk '{center=int(($2+$3)/2); print "chr"$1"\t"center"\t"center+1"\t"bed"`"NR"_dummy%"name}' bed=$3 name=$4 $1 > $2 
else
    awk '{center=int(($2+$3)/2); print $1"\t"center"\t"center+1"\t"bed"`"NR"_dummy%"name}' bed=$3 name=$4 $1 > $2 
    
fi 


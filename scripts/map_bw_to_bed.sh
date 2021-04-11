#!/bin/bash
# $1: input bed file
# $2: output csv gz file
# $3: slop 
# $4: bed file: just to get the slop
# $5: genoem file
# $6: kmeans cl id to select for
# $7: input bw_file 
bp_slop_bed=`echo $4 | awk -F'bp' '{print  $1}'`
bp_slop_bed=`echo 0 - $bp_slop_bed | bc`
bedtools slop -b ${bp_slop_bed} -g $5 -i $1 > ${1}.${6}.tmp.bed 
bedtools slop -b $3 -g $5 -i ${1}.${6}.tmp.bed | grep "\^${6}" > ${1}.${6}.${3}.tmp.bed 
python $NGS_SCRIPTS_DIR/map_bw_to_bed.py $7 ${1}.${6}.${3}.tmp.bed $2 

#!/bin/bash
# $1: bed file with motif center and unknown flank (will  determine here) 
# $2: bigwig file
# $3: output file
# $4: chip source key

awk '{print $1"\t"int(($2+$3)/2 - 300)"\t"int(($2+$3)/2 + 300)"\t"$4}' $1 > ${1}.${4}.tmp
python $NGS_SCRIPTS_DIR/map_bw_exact_stats_to_bed.py $2 ${1}.${4}.tmp $3


rm ${1}.${4}.tmp

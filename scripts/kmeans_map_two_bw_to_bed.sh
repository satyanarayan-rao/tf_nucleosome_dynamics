#!/bin/bash
# $1: bigwig_1
# $2: bigwig_2
# $3: bed_file
# $4: output gz file 1
# $5: output gz file 2

python $NGS_SCRIPTS_DIR/map_bw_to_bed.py $1 $3 $4
python $NGS_SCRIPTS_DIR/map_bw_to_bed.py $2 $3 $5

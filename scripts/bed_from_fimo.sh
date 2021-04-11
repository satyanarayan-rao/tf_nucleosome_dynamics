#!/bin/bash
grep -v "^#" $1 > ${1}.tmp 
python $NGS_SCRIPTS_DIR/fimo2bed.py --fimo ${1}.tmp --out $2

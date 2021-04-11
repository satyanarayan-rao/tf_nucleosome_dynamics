#!/bin/bash
# $1: fimo.tsv from fimo
# $2: fimo.bed

# skip header and remove all commments

head -1 ${1} > ${1}.header
awk 'NR>1' $1 | grep -v "^#"  > ${1}.tmp1 
cat ${1}.header ${1}.tmp1 > ${1}.tmp2
python $NGS_SCRIPTS_DIR/fimo2bed.py --fimo ${1}.tmp2 --out $2

# flank the file file 

# cleanup 
rm ${1}.tmp1  ${1}.tmp2 ${1}.header


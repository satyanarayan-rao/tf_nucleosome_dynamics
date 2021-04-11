#!/bin/bash
# $1: fimo tsv file
# $2: fimo bed file with top motif
# $3: annotation

head -1 $1 > ${1}.header
awk 'NR>1' $1 | grep -v "^#" | sort -k2,2 -k6,6gr > ${1}.sorted
python scripts/select_first_of_each.py ${1}.sorted > ${1}.top_for_each
cat ${1}.header ${1}.top_for_each > ${1}.top_with_header
python $NGS_SCRIPTS_DIR/fimo2bed_flexible.py --fimo ${1}.top_with_header --out ${1}.bed_to_center

# center on motif because different TFs can have diffferent width 
awk -F'[\t|]' '{center=int(($2 + $3)/2); print $1"\t"center"\t"center+1"\t"$4"|noMotif`"NR"_dummy%"annotation"_"$NF"\t"$6"\t"$7}' annotation=$3 ${1}.bed_to_center > ${2}

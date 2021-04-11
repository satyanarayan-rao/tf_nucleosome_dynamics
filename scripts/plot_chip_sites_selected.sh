#!/bin/bash
# $1: chip_tsv 
# $2: output file 
# $3: system name
# $4: pseudo chip value
# $5: cfDNA directory
# $6: chip source
# $7: reference for ks test 
# $8: bxplt_min_value
# $9: bxplt_max_value
# ${10}: dotted line color
# ${11}: ref box color
# ${12}: gnuplot color file - used for coloring the line; same will be used here for boxplots
# ${13}: source annotation
# ${14}: boxplot pdf file
# ${15}: figure title
# ${16}: select_site_string

echo ${17}
echo $7 1>&2
echo -e "Cluster\tMean.ChIP" > ${1}_${7}_${16}.tmp.header.tsv 
grep "${16}" $1 | awk -F'^' '{print $NF}' | cat ${1}_${7}_${16}.tmp.header.tsv - > ${1}_${7}_${16}.tmp.tsv
Rscript scripts/boxplot_chip.R ${1}_${7}_${16}.tmp.tsv 1 2 $2 $4 "${3}-${5}" "Cluster "  "Cluster" "Mean ChIP (${13})"  "${7}" $8 $9 ${10} ${11} ${12}  ${14} ${15}  ${17}

# cleanup
rm ${1}_${7}_${16}.tmp.header.tsv ${1}_${7}_${16}.tmp.tsv

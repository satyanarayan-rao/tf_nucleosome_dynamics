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
# ${16}: with baplot png
# ${17}: median_ordered tsv file
# ${18}: mean_tsv file - to do expected length calculation and order

echo $7 1>&2
echo -e "Cluster\tMean.ChIP" > ${1}_${7}_tmp.header.tsv 
awk -F'^' '{print $NF}' $1 | cat ${1}_${7}_tmp.header.tsv - > ${1}_${7}_tmp.tsv
#echo "Rscript scripts/boxplot_chip.R ${1}_${7}_tmp.tsv 1 2 $2 $4 "${3}-${5}" "Cluster "  "Cluster" "Mean ChIP \(${13}\)"  "${7}" $8 $9 ${10} ${11} ${12}  ${14} ${15} $1" 
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Command %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo "Rscript scripts/boxplot_chip.R ${1}_${7}_tmp.tsv 1 2 $2 $4 "${3}-${5}" "Cluster "  "Cluster" "Mean ChIP \(${13}\)"  "${7}" $8 $9 ${10} ${11} ${12}  ${14} ${15} $1 ${16} ${17} ${18} ${19} ${20}"  
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Command %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
Rscript scripts/boxplot_chip_expeted_len.R ${1}_${7}_tmp.tsv 1 2 $2 $4 "${3}-${5}" "Cluster "  "Cluster" "Mean ChIP (${13})"  "${7}" $8 $9 ${10} ${11} ${12}  ${14} ${15} $1 ${16} ${17} ${18} ${19} ${20}

# cleanup
rm ${1}_${7}_tmp.header.tsv ${1}_${7}_tmp.tsv

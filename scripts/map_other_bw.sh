#!/bin/bash
# $1: input bed tsv file
# $2: bigwig file
# $3: output raw csv file
# $4: slop
# $5: chip_source
# $6: genome size file
# $7: next bigwig annotation (mnase for example)
# first prepare the bed file from tsv, put the cluster id in name

#awk '{print $1}' $1 > ${1}.${5}.${7}.tmp_cl_id
#cut -f1 --complement $1 | paste -d'^' - ${1}.${5}.${7}.tmp_cl_id > ${1}.${5}.${7}.bed
#bedtools slop -b $4 -g ${6} -i ${1}.${5}.${7}.bed > ${1}-${5}-${7}-${4}.bed 
#
#python $NGS_SCRIPTS_DIR/map_bw_to_bed.py $2 ${1}-${5}-${7}-${4}.bed $3
## cleanup 
#rm ${1}.${5}.${7}.tmp_cl_id ${1}.${5}.${7}.bed 

# just reducing filename length
awk '{print $1}' $1 > ${1%_reassigned_bed.tsv}.${7}.tmp
cut -f1 --complement $1 | paste -d'^' - ${1%_reassigned_bed.tsv}.${7}.tmp > ${1%_reassigned_bed.tsv}.${7}.bed
bedtools slop -b $4 -g ${6} -i ${1%_reassigned_bed.tsv}.${7}.bed > ${1%_reassigned_bed.tsv}-${7}-${4}.bed 

python $NGS_SCRIPTS_DIR/map_bw_to_bed.py $2 ${1%_reassigned_bed.tsv}-${7}-${4}.bed $3
# cleanup 
rm ${1%_reassigned_bed.tsv}.${7}.tmp ${1%_reassigned_bed.tsv}.${7}.bed  

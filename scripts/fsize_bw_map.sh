#!/bin/bash
# $1: input bed tsv file
# $2: bigwig file
# $3: output raw csv file
# $4: slop
# $5: fragment_size
# $6: genome size file
# first prepare the bed file from tsv, put the cluster id in name

awk '{print $1}' $1 > ${1}.${5}.tmp_cl_id
cut -f1 --complement $1 | paste -d'^' - ${1}.${5}.tmp_cl_id > ${1}.${5}.bed 
bedtools slop -b $4 -g ${6} -i ${1}.${5}.bed > ${1}-${5}-${4}.bed  

python $NGS_SCRIPTS_DIR/map_bw_to_bed.py $2 ${1}-${5}-${4}.bed $3 

# cleanup
rm ${1}.${5}.tmp_cl_id ${1}.${5}.bed 

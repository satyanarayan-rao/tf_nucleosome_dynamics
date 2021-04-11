#!/bin/bash
# $1: input bed tsv file
# $2: bigwig file
# $3: output raw csv file
# $4: slop
# $5: fragment_size
# $6: manual cluster id location
# $7: reassinged clusted id bed tsv file
# $8: chip label
# $9: genome size file
# first prepare the bed file from tsv, put the cluster id in name


######## work on starting tsv file to re-assign cluster id - and then sort the file ##### 
python scripts/reassign_cl_id.py $1 $6 ${1}-${5}-${4}-${8}-reassigned.tsv
sort -k1,1nr ${1}-${5}-${4}-${8}-reassigned.tsv > $7
#########################################################################################
awk '{print $1}' $7 > ${1}.${5}.${4}.${8}.tmp_cl_id
cut -f1 --complement $7 | paste -d'^' - ${1}.${5}.${4}.${8}.tmp_cl_id > ${1}.${5}.${4}.${8}.bed 
bedtools slop -b $4 -g ${9} -i ${1}.${5}.${4}.${8}.bed >  ${1}.${5}.${4}.${8}-reassigned.bed
python $NGS_SCRIPTS_DIR/map_bw_to_bed.py $2 ${1}.${5}.${4}.${8}-reassigned.bed $3 

# cleanup
rm ${1}.${5}.${4}.${8}.tmp_cl_id ${1}-${5}-${4}-${8}-reassigned.tsv  

#!/bin/bash
# $1: tsv file contains entries where values were not found for `c` class fragments : assigned as cluster id 0
# $2: tsv file with cluster order 
# $3: motif_centered_bed
# $4: hg38.genome
# $5: slop val
# $6: lag val #

# FOR NOW: slop and lag val are same and chrom loci already slopped are required
#awk -F'[|\t]' '{print $2"`0"}' $1 > ${1}.tmp # assigning `0` for zeros file
#awk -F'[-:`]' '{print $1"\t"$2"\t"$3"\t"$0}' ${1}.tmp >  ${1}.tmp.bed
#bedtools slop -b $5 -g $4 -i ${1}.tmp.bed > ${1}.tmp.slopped.bed 

##### ---------------------- ##### 

#awk -F'[|\t]' '{print $2"`"$NF}' $2 > ${2}.tmp 
#awk -F'[-:`]' '{print $1"\t"$2"\t"$3"\t"$0}' ${2}.tmp >  ${2}.tmp.bed
#bedtools slop -b $5 -g $4 -i ${2}.tmp.bed > ${2}.tmp.slopped.bed
awk -F'[@\t]' '{print $1"\t"$2"\t"$3"\t"$4"`"$NF}' $2 > $3 

#cat ${1}.tmp.slopped.bed ${2}.tmp.slopped.bed > $3
#cat ${2}.tmp.slopped.bed > $3

### clean up

#rm ${1}.tmp ${1}.tmp.bed ${1}.tmp.slopped.bed  ${2}.tmp ${2}.tmp.bed ${2}.tmp.slopped.bed


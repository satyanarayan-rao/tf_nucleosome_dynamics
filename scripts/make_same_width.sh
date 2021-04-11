#!/bin/bash
# $1: wMotif file (slopped)
# $2: woMotif (not slopped : $3 - $2 = 1)
# $3: output file
# $4: genome size file

# first get the width of reference

width=`head -1 $1 | awk '{print $3 - $2}'`
# check if width is odd or even
v=`echo "$width % 2" | bc`

slop=`echo "${width}/2" | bc -l`
# if odd simply flank both side by ${width}/2 
if [ $v -eq 1 ]
then
     bedtools slop -b $slop -g $4 -i $2 > $3 
else
     bedtools slop -b $slop -g $4 -i $2 > ${3}.tmp
     cat ${3}.tmp | cut -f1-2 > ${3}.first_two
     cat ${3}.tmp | cut -f3 | awk '{print $1 - 1}' > ${3}.three
     cat ${3}.tmp | cut -f4- > ${3}.rest
     paste ${3}.first_two ${3}.three ${3}.rest > $3 
     # cleanup 
     rm ${3}.tmp ${3}.first_two ${3}.three ${3}.rest
    
fi 

#!/bin/bash
# $1: reassigned tsv
# $2: reassigned bed
awk '{print $1}' $1 > ${1}.tmp_cl_id
cut -f1 --complement $1 | paste -d'^' - ${1}.tmp_cl_id > $2

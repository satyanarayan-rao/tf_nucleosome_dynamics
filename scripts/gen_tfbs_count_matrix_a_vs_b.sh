#!/bin/bash
# $1: reassigned cluster tsv -> ih02 
# $2: reassigned cluster tsv -> mcf7_merged
# $3: matrix output

python scripts/table_find_tfbs_overlap.py $1 $2 | cut -f2- | sed '$d' | sed 's:\t$::g'> $3


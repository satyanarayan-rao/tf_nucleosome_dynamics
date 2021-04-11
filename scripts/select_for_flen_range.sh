#!/bin/bash
# $1: input count matrix: count matrix has 4 first columns as bed entries to traceback; 5th column is count for flen = 0; 6th for flen = 1; 7th for flen = 2; ... 
# $2: output matrix
# $3: first k column to select for (important - contain information)
# $4: min_flen (included)
# $5: max_flen (included)
flen_range_start=`echo "$3 + $4 + 1" | bc`
flen_range_end=`echo "$3 + $5 + 1" | bc`
cut -f1-${3},${flen_range_start}-${flen_range_end} $1 > $2


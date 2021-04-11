#!/bin/bash
# $1: input count matrix file
# $2: output kde matrix file
# $3: kde_n
# $4: kde_bw
# $5: min flen
# $5: max flen
Rscript scripts/kde_flen_range.R $1 $2 $3 $4 $5 $6

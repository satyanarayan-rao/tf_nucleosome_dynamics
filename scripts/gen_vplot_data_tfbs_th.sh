#!/bin/bash
# $1: cluster tagged intersection file
# $2: output file
# $3: start bed file input_bed/50b_* : just to get column numbers
# $4: bedgz length column
# $5: the flen_count_matrix/*flen_in_motif.tsv 
# $6: min_count to consider
# $7: flank in which vplot is being calculated
# $8: actual flank in the bed file
# $9: minimum fragment len towards counting the threshold
# ${10}: maximum fragment len towards counting the threhold
# ${11}: sites passing threshold
# get number of columns in bed

bed_cols=`awk '{print NF}' $3 | sort | uniq -c | awk '{print $2}'`
total=`echo "$bed_cols + $4" | bc -l`
python scripts/gen_vplot_data_tfbs_th.py $1 $total $bed_cols $2 $5 $6 $7  $8 ${9} ${10} ${11}

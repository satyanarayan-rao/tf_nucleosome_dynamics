#!/bin/bash
# $1: input matrix
# $2: input gnuplot to append in the beginning
# $3: eps file
# $4: pdf file 
# $5: output gplt file
# $6: title
# $7: bed file; just to count the lines
# $8: xlim
cp $2 $5 
tot_lines=`wc -l $7 | awk '{print $1}'`
tot_cols=`zless $1 | head -1 | awk '{print NF - 1}'`
cat << EOT >> $5
set output "$3"
set title "$6 ( n = ${tot_lines})" 
set yra [1:${tot_lines}]
set xra [1: ${tot_cols}]
set xtics ("-$8" 1, "0" $8, "+$8" $tot_cols) 
plot "<zcat $1" matrix with image notitle
EOT

gnuplot $5
ps2pdf  $3 $4


#!/bin/bash
# $1: input matrix
# $2: input bed kmeans # just to count total numbers in clusters and `total` lines
# $3: input base gnuplt file
# $4: heatmap gnuplot file
# $5: output eps 
# $6: output pdf 
# $7: sample name
cp $3 $4

lines_to_draw=`awk '{print $1}' $2 | uniq -c | awk '{sum+=$1; print sum}' | sed '$d'`
total_lines=`wc -l $2 | awk '{print $1}'`
flank_bp=`echo $1 | awk -F'_intersect_' '{print $NF}'| awk -F'bp' '{print $1}'`
cat <<EOT >> $4
set title '$7 (TFBS +- ${flank_bp} bp)'
set output '$5'
set yrange [1:$total_lines]
set xrange [1:100]
set cbra [0:0.05]
set xlabel 'Fragment length (bp)'
set xtics ("35" 1, "110" 51, "250" 100)
EOT
cnt=1
for line_loc in $lines_to_draw 
do
cat <<EOT >> $4
set arrow $cnt from 1, $line_loc to 100, $line_loc nohead lc rgb "#000000" front
EOT
cnt=`expr $cnt + 1`
done 

echo "plot '< cut -f1 --complement $1' matrix with image notitle" >> $4

gnuplot $4 
ps2pdf $5 $6

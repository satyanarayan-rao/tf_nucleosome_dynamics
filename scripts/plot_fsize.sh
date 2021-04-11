#!/bin/bash
# $1: input enrihcment over mean file 
# $2: bed tsv file: to get line locations 
# $3: base gnuplot file
# $4: output plt file 
# $5: output eps 
# $6: output pdf
# $7: sample name
# $8: slop value
# $9: fragment_size 
# ${10}: colmean file - used to set the cbra

cp $3 $4
lines_to_draw=`awk '{print $1}' $2 | uniq -c | awk '{sum+=$1; print sum}' | sed '$d'` 
total_lines=`wc -l $2 | awk '{print $1}'`
flank_bp=`echo $1 | awk -F'_intersect_' '{print $NF}'| awk -F'bp' '{print $1}'`
total_columns=`zcat $1 | head -1 | awk '{print NF-1}'`
left_x_tic=1
center_x_tic=`echo "$total_columns/2"| bc`
lr_x_tic_label=`echo $8 + $flank_bp | bc` 
cbra_min=`awk '{print $2}' ${10} | sort -g | head -1`
cbra_max=`awk '{print $2}' ${10} | sort -gr | head -1`

cat <<EOT >> $4 
set title '$7 (TFBS +- ${flank_bp} bp; $9)' 
set output '$5'
set yrange [1:$total_lines] 
set xrange [1:$total_columns]
set cbra [${cbra_min}:${cbra_max}]
set xtics ("-$lr_x_tic_label" 1, "0" $center_x_tic, "+$lr_x_tic_label" $total_columns) 
EOT
cnt=1
for line_loc in $lines_to_draw
do
cat <<EOT >> $4
set arrow $cnt from 1, $line_loc to $total_columns, $line_loc nohead lc rgb "#000000" front
EOT
cnt=`expr $cnt + 1`
done 

echo "plot '<zcat $1' matrix with image notitle" >> $4 

gnuplot $4 
## ps2pdf has file path restriction of 255 characters: 
## so first change to destionation dir and then do ps2pdf - then change back to dir

cwd=`pwd`
file_dir=`echo $5 | rev | cut -d'/' -f2- | rev` 
fname=`echo $5 | rev |   cut -d'/' -f1 | rev` 
cd $file_dir

ps2pdf $fname ${fname%.eps}.pdf  
cd $cwd

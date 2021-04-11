#!/bin/bash
# $1: gnuplot base file
# $2: gz file 1 
# $3: gz file 2
# $4: gnuplt_out_file
# $5: output eps file
# $6: output pdf file
# $7: slop value
# $8: fragment range `a`
# $9: fragment range `c`
# $10: directory name 
# $11: system name
# ${12}: frag label a 
# ${13}: frag label c
# ${14}: row_order file 

total_lines_a=`python $NGS_SCRIPTS_DIR/count_gz_file_lines.py $2`
total_lines_c=`python $NGS_SCRIPTS_DIR/count_gz_file_lines.py $3`
twice_of_slop=`echo ${7}*2 | bc`
cp $1 $4
total_lines=`wc -l ${14} | awk '{print $1}'`
awk '{print $2}' ${14} | uniq -c | awk '{print $2"\t"$1}' > ${14}.both_cl_id_counts
python $NGS_SCRIPTS_DIR/generate_cumsum.py ${14}.both_cl_id_counts > ${14}.both_cl_id_start_end
sed '$d' ${14}.both_cl_id_start_end  > ${14}.both_for_lines
cat <<EOD >> $4
set output '$5'
set ylabel "Motifs in clusters"
set xlabel "Distance from motif center (bp)"
set xtics ("-$7" 0, "0" $7, "+$7" $twice_of_slop)
set multiplot layout 1,2 rowsfirst title 'Fragment enrichment in clusters ${11} (${10})'
set title '${12} (n = $total_lines_a)' 
set yrange [1:$total_lines_a]
set cbra [1:4]
EOD
while read cl_id start end 
do
cat <<EOT >> $4
set arrow  $cl_id from 0, $end to $twice_of_slop, $end nohead lc rgb "#000000" front 
EOT
done < ${14}.both_for_lines
cat <<EOJ >> $4
plot '<zcat ${2}' matrix with image notitle 

EOJ
cat <<EOT >> $4
set title '${13} (n = $total_lines_c)'
set yrange [1:$total_lines_c]
set cbra [0.6:1.4]
EOT
while read cl_id start end 
do
cat <<EOT >> $4
set arrow  $cl_id from 0, $end to $twice_of_slop, $end nohead lc rgb "#000000" front 
EOT
done < ${14}.both_for_lines
cat <<EOJ >> $4
plot '<zcat ${3}' matrix with image notitle
EOJ

gnuplot $4
ps2pdf $5 $6

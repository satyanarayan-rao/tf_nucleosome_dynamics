#!/bin/bash

# $1: file to plot
# $2: gnuplot base file
# $3: gnuplot file 
# $4: plot eps
# $5: plot pdf
# $6: min_dist
# $7: max_dist
# $8: row_order file; just to count the lines
# $9: lag.val 
# ${10}: a
# ${11}: c 
# ${12 }: ccr title 
cp $2 $3

total_lines=`wc -l $8 | awk '{print $1}'`
awk '{print $2}' $8 | uniq -c | awk '{print $2"\t"$1}' > ${8}.cl_id_counts
python $NGS_SCRIPTS_DIR/generate_cumsum.py ${8}.cl_id_counts > ${8}.cl_id_start_end
sed '$d' ${8}.cl_id_start_end  > ${8}.for_lines
twice_lag_val=`echo ${9}*2 | bc` 
cat <<EOD >> $3
set output '$4'
set ylabel "Footprints in range (${6}-${7})"
set yrange [1:$total_lines] 
set xtics ("-${9}" 0, "0" ${9}, "${9}" $twice_lag_val)
EOD

# Draw lines after end of every cluster except last
while read cl_id start end 
do
cat <<EOT >> $3
set arrow  $cl_id from 0, $end to $twice_lag_val, $end nohead lc rgb "#000000" front 
EOT
done < ${8}.for_lines
cat <<EOD >> $3
set title '${12} (n = $total_lines)'
plot '<zcat $1' matrix with image notitle 
EOD

gnuplot $3
ps2pdf $4 $5 


# cleanup
rm ${8}.for_lines ${8}.cl_id_counts ${8}.cl_id_start_end

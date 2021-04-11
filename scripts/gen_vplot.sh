#!/bin/bash
# $1: vplot data
# $2: vplot gplt
# $3: vplot pdf
# $4: flank
# $5: sample name
# $6: bed file
# $7: cl id

# $8: zoomed gplt
# $9: zoomed pdf 
# ${10}: cbra limit 
# ${11}: reassigned_bed_tsv 
# ${12}: layered plot gplt
# ${13}: layered pdf 
# ${14}: zscore tsv
# ${15}: zscore gplt
# ${16}: zscore whole 
# ${17}: zscore split

min_frag_len=`head -1 $1 | awk '{print $1}'`
max_frag_len=`tail -n 1 $1 | awk '{print $1}'`
xra=`head -1 $1 | awk '{print NF -1}'`
eighty_on_plot=` echo "80 - $min_frag_len" | bc`
hundred_on_plot=` echo "100 - $min_frag_len" | bc`
onethirty_on_plot=` echo "130 - $min_frag_len" | bc`
onefifty_on_plot=` echo "150 - $min_frag_len" | bc`
oneeighty_on_plot=` echo "180 - $min_frag_len" | bc`
twohundred_on_plot=` echo "200 - $min_frag_len" | bc`
quarter_x_tic=`echo "${4}/2" |bc`
tfbs_half_width=`echo "($xra - 2*${4})/2" | bc`
nsites=`cat ${11} | awk '{if ($1==cl_id){print $0}}' cl_id=$7 | wc -l`

cat << EOF > $2
set terminal postscript enhanced color size 10, 5 font 'Arial, 15'
set size ratio 0.5
#set palette defined ( 0 "white", 0.25 "#c994c7", 0.375 "#78c679", 0.5 "#2c7fb8", 1 "#feb24c")
set palette defined ( 0 "#f2f0f7", 0.12 "#bdd7e7",  0.4 "#3182bd", 0.6 "#08519c", 0.8 "#9e9ac8", 1 "#756bb1")
#set cbra [0:${10}]
set output "${3%.pdf}.eps" 
set title "${5} on ${6} (cl${7}; n = ${nsites})"   noenhanced
set xra [0:$xra]
set yra[0:${twohundred_on_plot}]
#set ytics ("${min_frag_len}" 0, "80" ${eighty_on_plot}, "100" ${hundred_on_plot}, "130" ${onethirty_on_plot}, "150" ${onefifty_on_plot}, "180" ${oneeighty_on_plot}, "200" ${twohundred_on_plot}) 

set xtics ( "-${4}" 0,  "-${quarter_x_tic}" $quarter_x_tic, "0" `echo "${4} + ${tfbs_half_width}" | bc ` , "$quarter_x_tic" `echo "$4 + ${quarter_x_tic}"|bc`, "${4}" `echo "2*${4}" | bc` )  

plot '< cut -f2- $1' matrix with image notitle 
EOF

gnuplot $2

cwd=`pwd`
dest_dir=`echo $3 | rev | cut -d'/' -f2- | rev`
cd $dest_dir
fname=`echo $3 | rev | cut -d'/' -f1 | rev`

echo "ps2pdf ${fname%.pdf}.eps ${fname%.pdf}_tmp.pdf"
ps2pdf ${fname%.pdf}.eps ${fname%.pdf}_tmp.pdf

echo "pdfcrop ${3%.pdf}_tmp.pdf $3"
pdfcrop ${fname%.pdf}_tmp.pdf $fname

rm ${fname%.pdf}.eps ${fname%.pdf}_tmp.pdf
cd $cwd


########### Do zoomed v-pltos  ####### 
### -100 to 100  from TBSS center #### 
center_minus_hundred=`echo "$xra/2 - 100" | bc`
center_plus_hundred=`echo "$xra/2 + 100" | bc`


cat << EOF > $8
set terminal postscript enhanced color size 10, 5 font 'Arial, 15'
set size ratio 1
#set palette defined ( 0 "white", 0.25 "#c994c7", 0.375 "#78c679", 0.5 "#2c7fb8", 1 "#feb24c")
set palette defined ( 0 "#f2f0f7", 0.12 "#bdd7e7",  0.4 "#3182bd", 0.6 "#08519c", 0.8 "#9e9ac8", 1 "#756bb1")
#set cbra [0:${10}]
set output "${9%.pdf}_no_line.eps" 
set title "${5} on ${6} (cl${7}; n = ${nsites})"  noenhanced
set xra [${center_minus_hundred}:${center_plus_hundred}]
set yra[0:${twohundred_on_plot}]
#set ytics ("${min_frag_len}" 0, "80" ${eighty_on_plot}, "100" ${hundred_on_plot}, "130" ${onethirty_on_plot}, "150" ${onefifty_on_plot}, "180" ${oneeighty_on_plot}, "200" ${twohundred_on_plot}) 

set xtics ( "-100" ${center_minus_hundred},  "-50" `echo "${center_minus_hundred} + 50 " | bc -l`, "0" `echo "${xra}/2" | bc -l` , "+50" `echo "${xra}/2 + 50"|bc -l`, "+100" ${center_plus_hundred} ) 
set style line 1 lw 2 lc rgb 'green' 
set style line 2 lw 2 lc rgb 'red' 

plot '< cut -f2- $1' matrix with image notitle # without straight line y = 2*X

set output "${9%.pdf}_line.eps" 
plot '< cut -f2- $1' matrix with image notitle, 2*(x - `echo "($center_minus_hundred + $center_plus_hundred)/2" | bc -l`) + 10 ls 1 notitle, -2*(x -`echo "($center_minus_hundred + $center_plus_hundred)/2 " | bc -l`) + 10 ls 1 notitle, 2*(x - `echo "($center_minus_hundred + $center_plus_hundred)/2" | bc -l`) - 10 ls 2 notitle, -2*(x - `echo "($center_minus_hundred + $center_plus_hundred)/2 " | bc -l`) - 10 ls 2 notitle # with line Y = 2*X

EOF

gnuplot $8

cwd=`pwd`
dest_dir=`echo $9 | rev | cut -d'/' -f2- | rev`
cd $dest_dir
fname=`echo $9 | rev | cut -d'/' -f1 | rev`

echo "ps2pdf ${fname%.pdf}_no_line.eps ${fname%.pdf}_tmp.pdf"
ps2pdf ${fname%.pdf}_no_line.eps ${fname%.pdf}_tmp.pdf

echo "pdfcrop ${9%.pdf}_tmp.pdf $3"
pdfcrop ${fname%.pdf}_tmp.pdf ${fname} 
convert -density 300 -background white -alpha remove -trim ${fname} ${fname%.pdf}.png

echo "ps2pdf ${fname%.pdf}_line.eps ${fname%.pdf}_line.pdf" 

ps2pdf ${fname%.pdf}_line.eps ${fname%.pdf}_line.pdf 


convert -density 150 -background white -alpha remove -delay 150 -loop 2  ${fname} ${fname%.pdf}_line.pdf  ${fname%.pdf}.gif

rm ${fname%.pdf}_no_line.eps ${fname%.pdf}_tmp.pdf ${fname%.pdf}_line.eps
cd $cwd


#############  zoomed with splits - normal and zoomed ####### 

cat << EOF > ${12}
set terminal pdf enhanced color size 6, 1.5 font 'Arial, 15'
set size ratio 0.25
#set palette defined ( 0 "white", 0.25 "#c994c7", 0.375 "#78c679", 0.5 "#2c7fb8", 1 "#feb24c")
set palette defined ( 0 "#f2f0f7", 0.12 "#bdd7e7",  0.4 "#3182bd", 0.6 "#08519c", 0.8 "#9e9ac8", 1 "#756bb1")
#set cbra [0:${10}]
set cbtics format '%03.1f' left font 'Arial, 8'
set output "${13%.pdf}_0_50.pdf" 
#set title "${5} on ${6} (cl${7}; n = ${nsites})"  noenhanced
set xra [${center_minus_hundred}:${center_plus_hundred}]
set yra[0:50]
set ytics 25
set ytics format '%03.0f'
set noxtics

plot '< cut -f2- $1' matrix with image notitle 

set output "${13%.pdf}_51_100.pdf" 
set yra [51:100]
plot '< cut -f2- $1' matrix with image notitle 

set output "${13%.pdf}_101_150.pdf"  
set yra [101:150]
plot '< cut -f2- $1' matrix with image notitle 

set output "${13%.pdf}_151_200.pdf"  
set yra [151:200]
plot '< cut -f2- $1' matrix with image notitle 

EOF

gnuplot ${12}
pdfjam --trim '0cm 0cm 0.75cm 0.75cm' ${13%.pdf}_151_200.pdf  ${13%.pdf}_101_150.pdf  ${13%.pdf}_51_100.pdf ${13%.pdf}_0_50.pdf --nup 1x4 -o ${13} 
rm ${13%.pdf}_151_200.pdf  ${13%.pdf}_101_150.pdf  ${13%.pdf}_51_100.pdf ${13%.pdf}_0_50.pdf 

convert -density 300 -background white -alpha remove -trim ${13} ${13%.pdf}.tmp.png
convert ${13%.pdf}.tmp.png -background white -pointsize 75 label:"${5} on ${6} (cl${7}; n = ${nsites})"  +swap -gravity Center -append ${13%.pdf}.png

rm ${13%.pdf}.tmp.png


######### Do Z-scpre plotting ########### 
# first convert zero padded tsv to ######
cat $1 | cut -f2- > ${1}_onlydat.tsv
Rscript scripts/norm_tsv_matix.R ${1}_onlydat.tsv ${14}
rm ${1}_onlydat.tsv

cat <<EOF > ${15}
set terminal pdf enhanced color size 6, 6 font 'Arial, 15'
set size ratio 1
#set palette defined ( 0 "white", 0.25 "#c994c7", 0.375 "#78c679", 0.5 "#2c7fb8", 1 "#feb24c")
set palette defined ( 0 "#f2f0f7", 0.12 "#bdd7e7",  0.4 "#3182bd", 0.6 "#08519c", 0.8 "#9e9ac8", 1 "#756bb1")
#set cbra [0:${10}]
set cbtics format '%03.1f' left font 'Arial, 8'
set output "${16}" 
set title "${5} on ${6} (cl${7}; n = ${nsites})"  noenhanced
set xra [${center_minus_hundred}:${center_plus_hundred}]
set yra[0:200]
set ytics 25
set ytics format '%03.0f'
set xtics ( "-100" ${center_minus_hundred},  "-50" `echo "${center_minus_hundred} + 50 " | bc -l`, "0" `echo "${xra}/2" | bc -l` , "+50" `echo "${xra}/2 + 50"|bc -l`, "+100" ${center_plus_hundred} ) 

plot '${14}' matrix with image notitle 

set terminal pdf enhanced color size 6, 1.5 font 'Arial, 15'
set size ratio 0.25
set cbtics format '%03.1f' left font 'Arial, 8'
set output "${17%.pdf}_0_50.pdf" 
set yra [0:50]
set ytics 25
set ytics format '%03.0f'
set notitle
set noxtics

plot '${14}' matrix with image notitle 


set output "${17%.pdf}_51_100.pdf" 
set yra [51:100]
set ytics 25
set ytics format '%03.0f'
set noxtics
set notitle
plot '${14}' matrix with image notitle 

set output "${17%.pdf}_101_150.pdf" 
set yra [101:150]
set ytics 25
set ytics format '%03.0f'
set noxtics
set notitle
plot '${14}' matrix with image notitle 

set output "${17%.pdf}_151_200.pdf" 
set yra [151:200]
set ytics 25
set ytics format '%03.0f'
set notitle
set noxtics

plot '${14}' matrix with image notitle 

EOF

gnuplot ${15} 


convert -density 300 -background white -alpha remove -trim ${16} ${16%.pdf}.png
pdfjam --trim '0cm 0cm 0.75cm 0.75cm' ${17%.pdf}_151_200.pdf  ${17%.pdf}_101_150.pdf  ${17%.pdf}_51_100.pdf ${17%.pdf}_0_50.pdf --nup 1x4 -o ${17}
rm ${17%.pdf}_151_200.pdf  ${17%.pdf}_101_150.pdf  ${17%.pdf}_51_100.pdf ${17%.pdf}_0_50.pdf 

convert -density 300 -background white -alpha remove ${17} -trim ${17%.pdf}.tmp.png

convert ${17%.pdf}.tmp.png -background white -pointsize 75 label:"${5} on ${6} (cl${7}; n = ${nsites})"  +swap -gravity Center -append ${17%.pdf}.png

rm ${17%.pdf}.tmp.png

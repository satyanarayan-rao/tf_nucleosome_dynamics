
bed_file=`echo $1 | awk -F'/' '{print $NF}'`
new_bed=`echo $2 | awk -F'/' '{print $NF}'`
for i in ih02 mcf7_merged pt65_merged
do
    echo "zcat vplot_flen_count_matrix/${i}_intersect_${bed_file%.bed}_flen_in_motif.tsvverbos.tsv.gz | python scripts/subset_vplot_verbose.py $2 | gzip - > vplot_flen_count_matrix/${i}_intersect_${new_bed%.bed}_flen_in_motif.tsvverbos.tsv.gz"
done 

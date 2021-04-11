s=" "
for i in `cat metadata/tmp_sys.tsv`
do
    s=`cat manual_cluster_assignment_files/${i}_intersect_50bp_encode_accepted_chrom_chrY_removed_non_overlapping_dmatrix_flen_min_35_max_250_least_read_5_bw_3_nclust_6_max_iter_250_slop_300_chip_ENCFF324NQZ_rep2_ref_cluster_6_manual_order.expt_len.tsv | cut -f2 | awk 'NR>1' | paste -s -d '\t'`
    echo -e "$i\t$s"
done > tmp/pu1_expt_len.tsv

s=" "

for i in `cat metadata/tmp_sys.tsv`
do
    s=`cat manual_cluster_assignment_files/${i}_intersect_50bp_lyl1_non_overlapping_chrY_removed_dmatrix_flen_min_35_max_250_least_read_5_bw_3_nclust_6_max_iter_250_slop_300_chip_lyl1_gse63484_ref_cluster_6_manual_order.expt_len.tsv | cut -f2 | awk 'NR>1' | paste -s -d '\t'`
    echo -e "$i\t$s"
done > tmp/lyl1_expt_len.tsv



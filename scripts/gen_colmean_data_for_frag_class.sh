while read sys cls
do 
    len_cls=`echo $cls | tr ',' ' '`
    for j in $len_cls
    do 
        echo "colmean_savgol/${sys}_intersect_50bp_encode_accepted_chrom_chrY_removed_non_overlapping_dmatrix_flen_min_35_max_250_least_read_5_bw_3_nclust_6_max_iter_250_slop_1000_fsize_cf_30_${j}_chip_GM_ENCFF289XSX_ref_cluster_6_savgol_51_order_3_colmean_eom.tsv" 
    done 
done < metadata/breakdown_map.tsv

#!/bin/bash
# snakemake --dag -np --snakefile discover_foxa1_binding_patterns.smk plots/chip_box_plots/pt4_merged_intersect_25bp_unique_encode_and_carroll_with_strand_chrom_filtered_dmatrix_flen_min_35_max_250_least_read_5_nclust_4_max_iter_250_slop_300_chip_mcf7_ENCFF008ABX_bxplt.png --configfile configs/config.yaml | dot -Tsvg > trials_for_plots/qpg.svg
bed_file="unique_carroll_encode_and_kabos_with_strand_chrom_filtered"
for i in `cat metadata/systems.tsv | grep -v "^#" | sed '/^$/d'`
do
    for j in `cat metadata/tfbs_flanks.tsv | grep -v "^#" | sed '/^$/d'`
    do
        for k in `cat metadata/chip_sources.tsv | grep -v "^#" | sed '/^$/d'`
        do
            echo -e plots/chip_box_plots/${i}_intersect_${j}bp_${bed_file}_dmatrix_flen_min_35_max_250_least_read_5_nclust_4_max_iter_250_slop_300_chip_${k}_bxplt.png 
        done
    done
done

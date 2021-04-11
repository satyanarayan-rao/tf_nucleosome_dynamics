import os
import sys
include: "snakemakes/overall_length_distribution.smk"
include: "snakemakes/cfdna_fragments_in_tfbs.smk"
include: "snakemakes/kde_and_kmeans.smk"
include: "snakemakes/heatmap_and_mean_kde_plots.smk"
include: "snakemakes/heatmap_of_fragments_in_clsuters.smk"
include: "snakemakes/chip_analysis.smk"
include: "snakemakes/generate_neg_ctl_sites.smk"
include: "snakemakes/vplot_analysis.smk"
include: "snakemakes/excess_ratio_analysis.smk"

       

#rule foxa1_bed_intersect_with_pairs: 
#    input:
#        
#    params:
#    output:
#    shell:
#
#rule count_read_centers_in_range: 
#    input:
#        selected_range = "selected_fragment_range/{sample}_intersect_{bed}_flen_min_{min_len}_max_{max_len}.tsv"
#    params:
#        ignore_first_k_columns = 4 # python index besed (0-based)
#    output:
#        count_file = "fragment_center_counts_in_tfbs/{sample}_intersect_{bed}_flen_min_{min_len}_max_{max_len}_count.tsv"
#    shell:
#        "sh scripts/sum_tsv.sh {input.selected_range} {params.ignore_first_k_columns} {output.count_file}"
#
#rule kde_on_flen_range:
#    input:
#        selected_range = "selected_fragment_range/{sample}_intersect_{bed}_flen_min_{min_len}_max_{max_len}.tsv"
#    params:
#        kde_n = 100, 
#        kde_bw = 5
#    output:
#        kde_matrix = "kde_matrix/{sample}_intersect_{bed}_flen_min_{min_len}_max_{max_len}_kde.tsv"
#    shell: 
#        "sh scripts/kde.sh {input.selected_range} {output.kde_matrix} {params.kde_n} {params.kde_bw} {wildcards.min_len} {wildcards.max_len}"
#rule select_flen_range_for_kde: 
#    input:
#        count_matrix = "flen_count_matrix/{sample}_intersect_{bed}_flen_count.tsv" 
#    params:
#        first_k_char_columns = 4
#    output:
#        selected_range = "selected_fragment_range/{sample}_intersect_{bed}_flen_min_{min_len}_max_{max_len}.tsv"
#    shell:
#        "sh scripts/select_for_flen_range.sh {input.count_matrix} {output.selected_range} {params.first_k_char_columns} {wildcards.min_len} {wildcards.max_len}"



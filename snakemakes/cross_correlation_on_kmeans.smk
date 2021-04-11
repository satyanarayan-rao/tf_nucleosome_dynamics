rule prepare_read_count_matrix_for_ccr:
    input:
        out_kmeans_file_bed_tsv = "kmeans_on_kde/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_bed.tsv"
    params:
    output:
        bed_file = "ccr_beds/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}.bed"  
    shell:
        "awk '{{print $2\"\t\"$3\"\t\"$4\"\t\"$5\"^\"$1}}' {input.out_kmeans_file_bed_tsv} > {output.bed_file}"
def get_bw_file_name_top (wildcards): 
    param_key ="_".join([wildcards.sample, wildcards.bed, wildcards.min_flen,
                         wildcards.max_flen, wildcards.min_reads,
                         wildcards.bw, wildcards.nclust, wildcards.max_iter,
                         wildcards.kcl])
    fname = config["ccr_params"][param_key] 
    return fname
rule map_bw_to_slopped_bed_top: 
    input:
        bed_file = "ccr_beds/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}.bed",  
        bw_file = lambda wildcards: get_bw_file_name_top (wildcards) , 
        genome = "metadata/hg38.genome"
    params:
    output:
        enriched_bw_map = "ccr_bw_map/top_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}.csv.gz" 
    shell:
        "sh scripts/map_bw_to_bed.sh {input.bed_file} {output.enriched_bw_map}"
        " {wildcards.slop} {wildcards.bed} {input.genome} {wildcards.kcl} {input.bw_file}"  

rule map_bw_to_slopped_bed_nuc: 
    input:
        bed_file = "ccr_beds/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}.bed",  
        bw_file = lambda wildcards: config["flen_bigwigs"][wildcards.sample][wildcards.sample + "_125_170"],
        genome = "metadata/hg38.genome"
    params:
    output:
        enriched_bw_map = "ccr_bw_map/nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}.csv.gz" 
    shell:
        "sh scripts/map_bw_to_bed.sh {input.bed_file} {output.enriched_bw_map}"
        " {wildcards.slop} {wildcards.bed} {input.genome} {wildcards.kcl} {input.bw_file}"

rule crr_top_and_nuc: 
    input:
        enriched_bw_map_top = "ccr_bw_map/top_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}.csv.gz", 
        enriched_bw_map_nuc = "ccr_bw_map/nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}.csv.gz" 
    params:
    output: 
        top_ccr_nuc = "ccr_bw_map/ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}.csv.gz", 
        nan_out = "ccr_bw_map/nan_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}.tsv"
    shell:
        "Rscript scripts/cross_correlation_a_and_b.R {input.enriched_bw_map_top}"
        " {input.enriched_bw_map_nuc} {wildcards.lag}"
        " {output.top_ccr_nuc} {output.nan_out}"

rule kmeans_on_ccr: 
    input:
        top_ccr_nuc = "ccr_bw_map/ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}.csv.gz"
    params: 
    output: 
        kmeans_out = "kmeans_on_cluster_ccr/kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}.csv.gz", 
        row_order = "kmeans_on_cluster_ccr/kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}_row_order.tsv", 
    shell:
        "Rscript scripts/kmeans_cluster_eq_n.R {input.top_ccr_nuc}"
        " {output.kmeans_out} {output.row_order} {wildcards.kclust}"


rule plot_kmeans_on_ccr:
    input:
        kmeans_out = "kmeans_on_cluster_ccr/kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}.csv.gz", 
        row_order = "kmeans_on_cluster_ccr/kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}_row_order.tsv", 
        gnuplt_base_file = "utils/gnuplot_base_files/kmeans_dist.gplt"
    params: 
        ccr_label = "peak+-10bp vs nuc" 
    output: 
        kmeans_gnuplt_file = "plots/kmeans_on_cluster_ccr/kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}.gplt", 
        kmeans_plot_eps = "plots/kmeans_on_cluster_ccr/kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}.eps",
        kmeans_plot_pdf = "plots/kmeans_on_cluster_ccr/kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}.pdf",
    shell:
        "sh scripts/plot_kmeans.sh {input.kmeans_out} {input.gnuplt_base_file} {output.kmeans_gnuplt_file} {output.kmeans_plot_eps} {output.kmeans_plot_pdf} {wildcards.slop} {wildcards.slop} {input.row_order} {wildcards.slop} \"Peak\" \"Nuc\" \"{params.ccr_label}\""

rule generate_motif_centered_bed_with_slop: 
    input:
        nan_out = "ccr_bw_map/nan_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}.tsv", 
        row_order = "kmeans_on_cluster_ccr/kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}_row_order.tsv", 
        hg38_genome = "metadata/hg38.genome"
    output:
        motif_centered_bed = "beds_orderered_by_kmeans/kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}.bed"
    shell: 
        "sh scripts/generate_motif_centered_bed.sh {input.nan_out} {input.row_order}" 
        " {output.motif_centered_bed} {input.hg38_genome} {wildcards.slop} {wildcards.lag}"

rule generate_bw_kmeans_on_motif_centered_slop:
    input:
        motif_centered_bed = "beds_orderered_by_kmeans/kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}.bed", 
        peak_bw  = lambda wildcards: get_bw_file_name_top (wildcards), 
        nuc_bw = lambda wildcards: config["flen_bigwigs"][wildcards.sample][wildcards.sample + "_125_170"]
    params: 
    output:
        peak_gz = "bw_map_to_kmeans/kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}_peak.csv.gz", 
        peak_e_gz = "bw_map_to_kmeans/kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}_peak.csv.gz_e.csv.gz", 
        peak_exact_stats = "bw_map_to_kmeans/kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}_peak.csv.gz_exact_stats.tsv", 
        nuc_gz = "bw_map_to_kmeans/kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}_nuc.csv.gz", 
        nuc_e_gz = "bw_map_to_kmeans/kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}_nuc.csv.gz_e.csv.gz", 
        nuc_exact_stats = "bw_map_to_kmeans/kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}_nuc.csv.gz_exact_stats.tsv", 
       
    shell:
        "sh scripts/kmeans_map_two_bw_to_bed.sh {input.peak_bw} {input.nuc_bw}"
        " {input.motif_centered_bed} {output.peak_gz} {output.nuc_gz}"

rule plot_fragments_kmeans_on_ccr: 
    input:
        peak_gz = "bw_map_to_kmeans/kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}_peak.csv.gz_e.csv.gz", 
        nuc_gz = "bw_map_to_kmeans/kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}_nuc.csv.gz_e.csv.gz", 
        gnuplt_base_file = "utils/gnuplot_base_files/kmeans_frag_base.plt",
        row_order = "kmeans_on_cluster_ccr/kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}_row_order.tsv", 
        kmeans_plot_pdf = "plots/kmeans_on_cluster_ccr/kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}.pdf",
    params:
        system = "FoxA1",
        frag_label_peak = "35-345",
        frag_label_nuc = "125-170"
    output:
        kmeans_gnuplt_file = "plots/bw_map_to_kmeans/both_kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}.gplt",
        kmeans_plot_eps = "plots/bw_map_to_kmeans/both_kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}.eps", 
        kmeans_plot_pdf = "plots/bw_map_to_kmeans/both_kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}.pdf", 

    shell:
        "sh scripts/plot_kmeans_fragments.sh {input.gnuplt_base_file} {input.peak_gz} {input.nuc_gz} {output.kmeans_gnuplt_file} {output.kmeans_plot_eps} {output.kmeans_plot_pdf} {wildcards.slop} \"35-45\" \"125-170\" {wildcards.sample} {params.system} \"{params.frag_label_peak}\" \"{params.frag_label_nuc}\" {input.row_order}"

rule get_mean_chip_for_kmeans:
    input:
        motif_centered_bed = "beds_orderered_by_kmeans/kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}.bed", 
        chip_biwgwig = lambda wildcards: config["chip_bigwigs"][wildcards.chip]
    params:
    output:
        chip_exact_stats = "bw_map_to_kmeans/kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}_chip_{chip}.tsv",  
     
    shell: 
        "sh scripts/center_plus_mimus_300_chip.sh {input.motif_centered_bed} {input.chip_biwgwig} {output.chip_exact_stats} {wildcards.chip}" 


rule plot_mean_chip_for_kmeans: 
    input:
        chip_exact_stats = "bw_map_to_kmeans/kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}_chip_{chip}.tsv"
        
    params:
        chip_title = lambda wildcards: config["source_annotation"][wildcards.chip]
    output:
        chip_plot_ccr_within_cluster_png = "plots/bw_map_to_kmeans/chip_bxplt_kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}_chip_{chip}.png", 
        chip_plot_ccr_within_cluster_pdf = "plots/bw_map_to_kmeans/chip_bxplt_kmclust_{kclust}_ccr_top_vs_nuc_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slopped_{slop}_kmean_cl_{kcl}_lag_{lag}_chip_{chip}.pdf" 
    shell: 
        "Rscript scripts/plot_chip_within_cluster.R {input.chip_exact_stats} {output.chip_plot_ccr_within_cluster_png} {output.chip_plot_ccr_within_cluster_pdf} \"{params.chip_title}\""


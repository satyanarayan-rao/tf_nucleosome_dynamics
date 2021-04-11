rule flen_matrix_in_clsuters: 
    input:
        out_kmeans_file_bed_tsv = "kmeans_on_kde/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_nclust_{nclust}_max_iter_{max_iter}_bed.tsv", 
        bw_file = lambda wildcards: config["flen_bigwigs"][wildcards.sample][wildcards.sample + "_" + wildcards.fsize],
        genome_size_file = lambda wildcards: config["genome_size_map"][wildcards.sample]
    params:
    output:
        raw_signal = "fragment_bw_map/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}.csv.gz",
        enrichment_signal = "fragment_bw_map/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}.csv.gz_e.csv.gz",
        exact_signal = "fragment_bw_map/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}.csv.gz_exact_stats.tsv"
        
    shell:
        "sh scripts/fsize_bw_map.sh {input.out_kmeans_file_bed_tsv}"
        " {input.bw_file} {output.raw_signal} {wildcards.slop} {wildcards.fsize}"
        " {input.genome_size_file}"
rule plot_flen_matrix:
    input:
        enrichment_signal = "fragment_bw_map/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}.csv.gz_e.csv.gz", 
        out_kmeans_file_bed_tsv = "kmeans_on_kde/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_nclust_{nclust}_max_iter_{max_iter}_bed.tsv", # this will be used just for cluster lines on the heatmap 
        base_gnuplt_file = "utils/gnuplot_base_files/fragment_matrix.gplt"
        
    params:
        
    output:
        eom_eps = "plots/fragment_plots/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_eom.eps", 
        eom_pdf = "plots/fragment_plots/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_eom.pdf",
        eom_plt = "plots/fragment_plots/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_eom.gplt"
    shell: 
        "sh scripts/plot_fsize.sh {input.enrichment_signal} {input.out_kmeans_file_bed_tsv}" 
        " {input.base_gnuplt_file} {output.eom_plt} {output.eom_eps} {output.eom_pdf} {wildcards.sample} {wildcards.slop} {wildcards.fsize}" 

########################################################################

########################################################################


def get_manual_cl_id_file_heatmap(wildcards): 
    # manually assinged cluster id can be found in `manual_cluster_assignment_files/`
    # with file names as <sample>_<bed>_<min_flen>_<max_flen>_<min_reads>_<bw>_<nclust>_<max_iter>_cl_id.tsv
    fname = "_".join([wildcards.sample, wildcards.bed, wildcards.min_flen,
                      wildcards.max_flen, wildcards.min_reads, wildcards.bw,
                      wildcards.nclust, wildcards.max_iter])
    fname = "manual_cluster_assignment_files/" + fname + ".tsv"
    return fname


rule flen_matrix_in_clsuters_manual: 
    input:
        out_kmeans_file_bed_tsv = "kmeans_on_kde/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_bed.tsv", 
        bw_file = lambda wildcards: config["flen_bigwigs"][wildcards.sample][wildcards.sample + "_" + wildcards.fsize],
        #manual_cl_id_file = lambda wildcards: get_manual_cl_id_file_heatmap(wildcards)
        #manual_cl_id_file = "manual_cluster_assignment_files/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_300_chip_{chip}_ref_cluster_{ref_for_ks_test}_manual_order.tsv",  #first get the order by chip and then use if for the manual plot   
        manual_cl_id_file = lambda wildcards: "manual_cluster_assignment_files/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_" + str(config["chip_bxplt_tfbs_flank"][wildcards.chip]) + "_chip_{chip}_ref_cluster_{ref_for_ks_test}_manual_order.tsv",  #first get the order by chip and then use if for the manual plot  
        genome_size_file = lambda wildcards: config["genome_size_map"][wildcards.sample]
        
    params:
    output:
        raw_signal = "fragment_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}.csv.gz",
        enrichment_signal = "fragment_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}.csv.gz_e.csv.gz",
        exact_signal = "fragment_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}.csv.gz_exact_stats.tsv", 

        reassigned_bed_tsv = "fragment_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_reassigned_bed.tsv"

    shell:
        "sh scripts/fsize_bw_map_manual.sh {input.out_kmeans_file_bed_tsv}"
        " {input.bw_file} {output.raw_signal}"
        " {wildcards.slop} {wildcards.fsize} {input.manual_cl_id_file} {output.reassigned_bed_tsv}"
        " {wildcards.chip} {input.genome_size_file}"
def get_gnuplot_base_file (wildcards): 
    fname = "utils/gnuplot_base_files/fragment_matrix_{s}.gplt".format(s=wildcards.fsize)
    return fname
rule plot_flen_matrix_manual:
    input:
        enrichment_signal = "fragment_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}.csv.gz_e.csv.gz",
        reassigned_bed_tsv = "fragment_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_reassigned_bed.tsv",
        base_gnuplt_file = lambda wildcards: get_gnuplot_base_file(wildcards),
        colmean_eom = "colmean_fragment_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_colmean_eom.tsv",
        
    params:
        
    output:
        eom_eps = "plots/fragment_plots_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_eom.eps", 
        eom_pdf = "plots/fragment_plots_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_eom.pdf",
        eom_plt = "plots/fragment_plots_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_eom.gplt"
    shell: 
        "sh scripts/plot_fsize.sh {input.enrichment_signal} {input.reassigned_bed_tsv}" 
        " {input.base_gnuplt_file} {output.eom_plt} {output.eom_eps} {output.eom_pdf}"
        " {wildcards.sample} {wildcards.slop} {wildcards.fsize} {input.colmean_eom}" 

rule cluster_wise_colmeans_data:
    input:
        enrichment_signal = "fragment_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}.csv.gz_e.csv.gz",
        colmean_line_colors = "utils/gnuplot_base_files/boxplot_colors.yaml" 
        
    params:
        total_slop = lambda wildcards: str(int(wildcards.slop) + int(wildcards.bed.split("bp_")[0]))
    output:
        colmean_eom = "colmean_fragment_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_colmean_eom.tsv",
    shell:
        "sh scripts/calculate_colmeans_of_mnase_eom.sh {input.enrichment_signal}"
        " {wildcards.nclust} {params.total_slop} {output.colmean_eom}"
rule cluster_wise_colmeans_plot:
    input:
        colmean_eom = "colmean_fragment_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_colmean_eom.tsv",
        colmean_line_colors = "utils/gnuplot_base_files/boxplot_colors.yaml" 
        
    params:
        total_slop = lambda wildcards: str(int(wildcards.slop) + int(wildcards.bed.split("bp_")[0]))
    output:
        colmean_eom_pdf = "plots/colmean_fragment_plots_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_colmean_eom.pdf",
        colmean_eom_png = "plots/colmean_fragment_plots_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_colmean_eom.png",
    shell:
        "Rscript scripts/heatmap_colmean_plot.R {input.colmean_eom}" 
        " {output.colmean_eom_pdf} {output.colmean_eom_png} {input.colmean_line_colors} {wildcards.nclust}"

################################## Savgol Smoothing and plotting ######################### 

rule savgol_smoothing:
    input:
        enrichment_signal = "fragment_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}.csv.gz_e.csv.gz",
    output: 
        savgol_smooth = "savgol_smooth_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_savgol_{window}_order_{order}.csv.gz_e.csv.gz"
    shell:
        "python $NGS_SCRIPTS_DIR/savgol_smooth.py {input.enrichment_signal} {output.savgol_smooth} {wildcards.window} {wildcards.order}"

rule colmean_savgol:
    input:
        savgol_smooth = "savgol_smooth_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_savgol_{window}_order_{order}.csv.gz_e.csv.gz"
    params:
        total_slop = lambda wildcards: str(int(wildcards.slop) + int(wildcards.bed.split("bp_")[0]))
    output:
        colmean_eom = "colmean_savgol/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_savgol_{window}_order_{order}_colmean_eom.tsv"
    shell:
        "sh scripts/calculate_colmeans_of_mnase_eom.sh {input.savgol_smooth}"
        " {wildcards.nclust} {params.total_slop} {output.colmean_eom}" 

rule savgol_cluster_wise_colmeans_plot:
    input:
        colmean_eom = "colmean_savgol/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_savgol_{window}_order_{order}_colmean_eom.tsv",
        colmean_line_colors = "utils/gnuplot_base_files/boxplot_colors.yaml" 
        
    params:
        total_slop = lambda wildcards: str(int(wildcards.slop) + int(wildcards.bed.split("bp_")[0]))
    output:
        colmean_eom_pdf = "plots/colmean_savgol_plots_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_savgol_{window}_order_{order}_colmean_eom.pdf",
        colmean_eom_png = "plots/colmean_savgol_plots_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_savgol_{window}_order_{order}_colmean_eom.png",
    shell:
        "Rscript scripts/heatmap_colmean_plot.R {input.colmean_eom}" 
        " {output.colmean_eom_pdf} {output.colmean_eom_png} {input.colmean_line_colors} {wildcards.nclust}"

rule with_and_without_savgol_plots:
    input:
        without_savgol = "colmean_fragment_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_colmean_eom.tsv",
        with_savgol = "colmean_savgol/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_savgol_{window}_order_{order}_colmean_eom.tsv", 
        colmean_line_colors = "utils/gnuplot_base_files/boxplot_colors.yaml"  
    output:
        compare_pdf = "plots/with_and_without_savgol/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_savgol_{window}_order_{order}_vs_nosavgol.pdf",
        compare_png = "plots/with_and_without_savgol/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_savgol_{window}_order_{order}_vs_nosavgol.png" 
    shell:
        "Rscript scripts/compare_colmean_nosavgol_vs_savgol.R {input.without_savgol}"
        " {input.with_savgol} {output.compare_pdf} {output.compare_png} {input.colmean_line_colors} {wildcards.nclust} {wildcards.sample}" 



############### plotting savgol heatmaps ###########
rule savgol_heatmap:
    input:
        enrichment_signal = "savgol_smooth_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_savgol_{window}_order_{order}.csv.gz_e.csv.gz",
        reassigned_bed_tsv = "fragment_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_reassigned_bed.tsv",
        base_gnuplt_file = lambda wildcards: get_gnuplot_base_file(wildcards),
        colmean_eom = "colmean_savgol/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_savgol_{window}_order_{order}_colmean_eom.tsv",
        
    params:
        
    output:
        eom_eps = "plots/savgol_plots_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_savgol_{window}_order_{order}.eps", 
        eom_pdf = "plots/savgol_plots_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_savgol_{window}_order_{order}.pdf",
        eom_plt = "plots/savgol_plots_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_savgol_{window}_order_{order}.gplt"
    shell: 
        "sh scripts/plot_fsize.sh {input.enrichment_signal} {input.reassigned_bed_tsv}" 
        " {input.base_gnuplt_file} {output.eom_plt} {output.eom_eps} {output.eom_pdf}"
        " {wildcards.sample} {wildcards.slop} {wildcards.fsize} {input.colmean_eom}"  

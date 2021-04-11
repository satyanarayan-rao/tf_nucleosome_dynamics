def get_figure_title(wildcards):
    tfbs_flank = "".join(["TFBS", u"\u00B1", wildcards.bed.split("_")[0]])
    bw_str = "=".join(["bw", wildcards.bw])
    title = "-".join([config["system"],
                      config["sample_title"][wildcards.sample],
                      bw_str, tfbs_flank])
    return title
rule mean_plot_of_clusters: 
    input:
        out_kmeans_file = "kmeans_on_kde/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}.tsv", 
        x_tics = "metadata/x_tics.tsv",
        gnuplot_color_file = "utils/gnuplot_base_files/boxplot_colors.yaml"
        
    params:
        figure_title = lambda wildcards: get_figure_title(wildcards)
    output:
        mean_line_plot_of_clusters = "plots/mean_plots_of_clusters/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_line_plot.png",
        mean_line_plot_of_clusters_pdf = "plots/mean_plots_of_clusters/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_line_plot.pdf",
        mean_tsv = "kde_cluster_colmeans/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}.tsv"
    shell:
        "Rscript scripts/line_plots_kmeans.R {input.out_kmeans_file} {output.mean_line_plot_of_clusters} {input.x_tics} {output.mean_line_plot_of_clusters_pdf} {params.figure_title} {input.gnuplot_color_file} {output.mean_tsv}" 

############################################################
#### Manually assign cluster id by recording appropriate ###
#### cluster id in file with visual inspection #############
############################################################
def get_manual_cl_id_file(wildcards): 
    # manually assinged cluster id can be found in `manual_cluster_assignment_files/`
    # with file names as <sample>_<bed>_<min_flen>_<max_flen>_<min_reads>_<bw>_<nclust>_<max_iter>_cl_id.tsv
    fname = "_".join([wildcards.sample, wildcards.bed, wildcards.min_flen,
                      wildcards.max_flen, wildcards.min_reads, wildcards.bw,
                      wildcards.nclust, wildcards.max_iter])
    fname = "manual_cluster_assignment_files/" + fname + ".tsv"
    return fname
rule manual_mean_plots_of_clusters: 
    input:
        mean_tsv = "kde_cluster_colmeans/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}.tsv", 
        manual_cl_id_file = lambda wildcards: get_manual_cl_id_file (wildcards),
        gnuplot_color_file = "utils/gnuplot_base_files/boxplot_colors.yaml"
    params: 
        figure_title = lambda wildcards: get_figure_title(wildcards)
    output:
        mean_line_plot_of_manual_clusters = "plots/manually_assigned_mean_plots_of_clusters/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_manual_line_plot.png", 
        mean_line_plot_of_manual_clusters_pdf = "plots/manually_assigned_mean_plots_of_clusters/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_manual_line_plot.pdf", 
        manual_mean_tsv = "kde_cluster_colmeans_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_manual.tsv",
        peak_location_tsv = "kde_cluster_colmeans_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_peak_location.tsv"
    shell:
        "Rscript scripts/line_plots_kmeans_manual.R {input.mean_tsv} {input.manual_cl_id_file} {output.mean_line_plot_of_manual_clusters} {params.figure_title} {input.gnuplot_color_file} {output.manual_mean_tsv} {output.peak_location_tsv} {output.mean_line_plot_of_manual_clusters_pdf}" 

def get_cluster_wise_files (wildcards):
    file_list = [] 
    for s in config["cluster_wise_mean_plots"][wildcards.setting]["systems"]:
        if s in ["mcf7_merged_mm10", "KH182_MsY"]:
            fname = "kde_cluster_colmeans/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}"\
                    "_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}"\
                    "_max_iter_{max_iter}.tsv".format(sample=s, 
                    bed="50bp_unique_carroll_encode_with_strand_chrom_filtered_mm10", 
                    min_flen=wildcards.min_flen, 
                    max_flen=wildcards.max_flen, min_reads = wildcards.min_reads,
                    bw=wildcards.bw, nclust = wildcards.nclust, 
                    max_iter = wildcards.max_iter)
        else:
            fname = "kde_cluster_colmeans/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}"\
                    "_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}"\
                    "_max_iter_{max_iter}.tsv".format(sample=s, 
                    bed=wildcards.bed, min_flen=wildcards.min_flen, 
                    max_flen=wildcards.max_flen, min_reads = wildcards.min_reads,
                    bw=wildcards.bw, nclust = wildcards.nclust, 
                    max_iter = wildcards.max_iter)
        file_list.append (fname)
    return file_list
rule cluster_wise_mean_plots:
    input: 
        input_files = lambda wildcards: get_cluster_wise_files(wildcards)
    params:
        cluster_id = lambda wildcards: config["cluster_wise_mean_plots"][wildcards.setting]["cluster_id"],
        input_files_label = lambda wildcards: config["cluster_wise_mean_plots"][wildcards.setting]["labels"],
        tf_system = config["system"]
    output:
        cluster_specific_plot = "plots/mean_plots_of_clusters/setting_{setting}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_grouped_plot.png"
    shell:
        "Rscript scripts/cluster_specific_flen_density_plot.R \"{input.input_files}\" {params.cluster_id} \"{params.input_files_label}\""
        " {output.cluster_specific_plot} {params.tf_system}"

rule heatmap_of_clusters: 
    input:
        out_kmeans_file = "kmeans_on_kde/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}.tsv", 
        out_kmeans_file_bed_tsv = "kmeans_on_kde/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_bed.tsv",
        base_gnuplt_file = "utils/gnuplot_base_files/kmeans_kde_heatmap.plt"
    params:
        system = config["system"]
    output:
        heatmap_file_eps = "plots/kmeans_kde_heatmap/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_heatmap.eps", 
        heatmap_file_pdf = "plots/kmeans_kde_heatmap/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_heatmap.pdf", 
        heatmap_file_plt = "plots/kmeans_kde_heatmap/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_heatmap.plt" 
    shell: 
        "sh scripts/kmeans_kde_heatmap.sh {input.out_kmeans_file} {input.out_kmeans_file_bed_tsv} {input.base_gnuplt_file} {output.heatmap_file_plt} {output.heatmap_file_eps} {output.heatmap_file_pdf} {wildcards.sample}"
 

import yaml
import re
rule calculate_chip_in_cluster:
    input:    
        out_kmeans_file_bed_tsv = "kmeans_on_kde/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_bed.tsv", 
        bw_file = lambda wildcards: config["chip_bigwigs"][wildcards.chip],
        genome_size_file = lambda wildcards: config["genome_size_map"][wildcards.sample]
    params:
    output:
        raw_signal = "bw_map/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_chip_{chip}.csv.gz",
        enrichment_signal = "bw_map/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_chip_{chip}.csv.gz_e.csv.gz",
        exact_signal = "bw_map/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_chip_{chip}.csv.gz_exact_stats.tsv"
    shell:
        "sh scripts/map_bw.sh {input.out_kmeans_file_bed_tsv} {input.bw_file} {output.raw_signal} {wildcards.slop} {wildcards.chip} {input.genome_size_file}"

def get_boxplot_title(wildcards):
    #tfbs_flank = "".join(["TFBS", u"\u00B1", wildcards.bed.split("_")[0]])
    #bw_str = "=".join(["bw", wildcards.bw])
    #sys_title = config["bed_annotation"][wildcards.bed[wildcards.bed.find("_")+1:len(wildcards.bed)]]
    #title = "-".join([ sys_title,#config["system"],
    #                  config["sample_title"][wildcards.sample],
    #                  bw_str, tfbs_flank])
    title = config["simplified_labels"][wildcards.bed[wildcards.bed.find("_")+1:len(wildcards.bed)]] + "-" + config["sample_title"][wildcards.sample]
    return title 
rule plot_chip_scores: 
    input:
        mean_line_plot_of_clusters = "plots/mean_plots_of_clusters/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_line_plot.png",
        heatmap_file_eps = "plots/kmeans_kde_heatmap/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_heatmap.eps", 
        exact_signal = "bw_map/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_chip_{chip}.csv.gz_exact_stats.tsv",
        gnuplot_color_file = "utils/gnuplot_base_files/boxplot_colors.yaml", 
        mean_tsv = "kde_cluster_colmeans/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}.tsv"

    params:
        #system = config["system"], 
        system = lambda wildcards: config["bed_annotation"][wildcards.bed[wildcards.bed.find("_")+1:len(wildcards.bed)]],
        pseudo_chip = 0.01,
	      ref_for_ks_test = lambda wildcards: wildcards.ref_for_ks_test,
        bxplt_min_value = lambda wildcards: config["plot_chip_params"][wildcards.bed]["min_value"],
        bxplt_max_value = lambda wildcards: config["plot_chip_params"][wildcards.bed]["max_value"],
        plt_width = lambda wildcards: config["plot_chip_params"][wildcards.bed]["width"],
        plt_height = lambda wildcards: config["plot_chip_params"][wildcards.bed]["height"],
        dotted_line_color = "\"#636363\"",
        ref_ks_box_color = lambda wildcards: config["cluster_colors"][wildcards.ref_for_ks_test], 
        source_annotation = lambda wildcards: config["source_annotation"][wildcards.chip], 
        sample_title = lambda wildcards: config["sample_title"][wildcards.sample],
        plot_title = lambda wildcards: get_boxplot_title(wildcards),
        #figure_title = lambda wildcards: get_figure_title(wildcards), 
        figure_title = lambda wildcards: "-".join([config["bed_annotation"][wildcards.bed[wildcards.bed.find("_")+1:len(wildcards.bed)]],
                          config["sample_title"][wildcards.sample], "bw=" + wildcards.bw, "".join(["TFBS", u"\u00B1", wildcards.bed.split("_")[0]]) ] ),
        log_flag = lambda wildcards: config["do_log"][wildcards.chip]
    output:
        chip_box_plot = "plots/chip_box_plots/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_chip_{chip}_ref_cluster_{ref_for_ks_test}_bxplt.png", 
        chip_box_plot_pdf = "plots/chip_box_plots/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_chip_{chip}_ref_cluster_{ref_for_ks_test}_bxplt.pdf",
        ratio_barplot = "plots/chip_box_plots/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_chip_{chip}_ref_cluster_{ref_for_ks_test}.with.barplot.png",
        chip_meadian_reorder_tsv = "manual_cluster_assignment_files/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_chip_{chip}_ref_cluster_{ref_for_ks_test}_manual_order.tsv",  #first get the order by chip and then use if for the manual plot  
        expt_len_file = "manual_cluster_assignment_files/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_chip_{chip}_ref_cluster_{ref_for_ks_test}_manual_order.expt_len.tsv",
        expt_verbose_file = "manual_cluster_assignment_files/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_chip_{chip}_ref_cluster_{ref_for_ks_test}_manual_order.expt_verbose.tsv", 
        

        # for mean line plots with manual clusters
        mean_line_plot_of_manual_clusters = "plots/manually_assigned_mean_plots_of_clusters/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_chip_{chip}_ref_cluster_{ref_for_ks_test}_manual_line_plot.png", 
        mean_line_plot_of_manual_clusters_pdf = "plots/manually_assigned_mean_plots_of_clusters/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_chip_{chip}_ref_cluster_{ref_for_ks_test}_manual_line_plot.pdf",
        manual_mean_tsv = "kde_cluster_colmeans_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_chip_{chip}_ref_cluster_{ref_for_ks_test}_manual.tsv",
        peak_location_tsv = "kde_cluster_colmeans_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_chip_{chip}_ref_cluster_{ref_for_ks_test}_peak_location.tsv",  
        median_ratio_plot = "plots/chip_box_plots/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_chip_{chip}_ref_cluster_{ref_for_ks_test}_median_ratio.pdf", 
        median_ratio_tsv = "plots/chip_box_plots/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_chip_{chip}_ref_cluster_{ref_for_ks_test}_median_ratio.tsv", 

    shell:
        #"Rscript scripts/boxplot_tsv.R {input.exact_signal} {output.chip_box_plot}"
        "sh scripts/plot_chip.sh {input.exact_signal} {output.chip_box_plot}"
        " {params.system} {params.pseudo_chip} {params.sample_title}"
        " {wildcards.chip} \"{params.ref_for_ks_test}\" \"{params.bxplt_min_value}\""
        " \"{params.bxplt_max_value}\" {params.dotted_line_color}"
        " {params.ref_ks_box_color} {input.gnuplot_color_file}"
        " \"{params.source_annotation}\" {output.chip_box_plot_pdf}"
        " \"{params.plot_title}\" {output.ratio_barplot} {output.chip_meadian_reorder_tsv}"
        " {input.mean_tsv} {output.expt_len_file} {output.expt_verbose_file} {params.log_flag}"
        " {output.median_ratio_plot} {params.plt_width} {params.plt_height} {output.median_ratio_tsv}"
        "; Rscript scripts/line_plots_kmeans_manual.R {input.mean_tsv}"
        " {output.chip_meadian_reorder_tsv} {output.mean_line_plot_of_manual_clusters}"
        " {params.figure_title} {input.gnuplot_color_file} {output.manual_mean_tsv}"
        " {output.peak_location_tsv} {output.mean_line_plot_of_manual_clusters_pdf}"
        " {output.expt_len_file}"



### Mapping MNase on "ChIP" ordered clusters ### 

### Currently clusters are ordered based on weighted lengths, but they still ## 
### go through ChIP based pipeline, so the workflow is like the following ### 
### ChIP ---> decides order -> bed files are oreded ---> heatmaps ### 
### Thus, MNase can't be treated as any other ChIP, but MNase heatmpas will depend on ChIP ### 

rule map_other_bigwigs:
    input:    
        reassigned_bed_tsv = "fragment_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_reassigned_bed.tsv", # taken from heatmap_of_fragments_in_clsuters.smk
        bw_file = lambda wildcards: config["chip_bigwigs"][wildcards.next_bw],
        genome_size_file = lambda wildcards: config["genome_size_map"][wildcards.sample]
    params:
    output:
        raw_signal = "other_bw_map/{sample}_intersect_{bed}_dm_fmin_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_{ref_for_ks_test}_oth_{next_bw}.csv.gz",
        enrichment_signal = "other_bw_map/{sample}_intersect_{bed}_dm_fmin_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_{ref_for_ks_test}_oth_{next_bw}.csv.gz_e.csv.gz",
        exact_signal = "other_bw_map/{sample}_intersect_{bed}_dm_fmin_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_{ref_for_ks_test}_oth_{next_bw}.csv.gz_exact_stats.tsv"
    shell:
        "sh scripts/map_other_bw.sh {input.reassigned_bed_tsv} {input.bw_file} {output.raw_signal} {wildcards.slop} {wildcards.chip} {input.genome_size_file} {wildcards.next_bw}"

rule cluster_wise_colmeans_data_other_bigwig:
    input:
        enrichment_signal = "other_bw_map/{sample}_intersect_{bed}_dm_fmin_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_{ref_for_ks_test}_oth_{next_bw}.csv.gz_e.csv.gz",
        #enrichment_signal = "other_bw_map/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_another_bw_{next_bw}.csv.gz_e.csv.gz",
        
    params:
        total_slop = lambda wildcards: str(int(wildcards.slop) + int(wildcards.bed.split("bp_")[0]))
    output:
        colmean_eom = "colmean_other_bw_map_manual/{sample}_intersect_{bed}_dm_fmin_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_{ref_for_ks_test}_oth_{next_bw}_colmean_eom.tsv",
    shell:
        "sh scripts/calculate_colmeans_of_mnase_eom.sh {input.enrichment_signal}"
        " {wildcards.nclust} {params.total_slop} {output.colmean_eom}"

def get_gnuplot_base_file_for_other_bw (wildcards): 
    fname = "utils/gnuplot_base_files/bigwig_{s}.gplt".format(s=wildcards.next_bw)
    return fname
rule plot_other_bigwig_manual:
    input:
        #enrichment_signal = "other_bw_map/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_another_bw_{next_bw}.csv.gz_e.csv.gz",
        enrichment_signal = "other_bw_map/{sample}_intersect_{bed}_dm_fmin_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_{ref_for_ks_test}_oth_{next_bw}.csv.gz_e.csv.gz",
        reassigned_bed_tsv = "fragment_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_reassigned_bed.tsv",
        base_gnuplt_file = lambda wildcards: get_gnuplot_base_file_for_other_bw(wildcards),
        #colmean_eom = "colmean_other_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_another_bw_{next_bw}_colmean_eom.tsv",
        colmean_eom = "colmean_other_bw_map_manual/{sample}_intersect_{bed}_dm_fmin_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_{ref_for_ks_test}_oth_{next_bw}_colmean_eom.tsv",
        
    params:
        
    output:
        eom_eps = "plots/other_bigwig_plots_manual/{sample}_intersect_{bed}_dm_fmin_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_{ref_for_ks_test}_oth_{next_bw}_eom.eps", 
        eom_pdf = "plots/other_bigwig_plots_manual/{sample}_intersect_{bed}_dm_fmin_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_{ref_for_ks_test}_oth_{next_bw}_eom.pdf",
        eom_plt = "plots/other_bigiwg_plots_manual/{sample}_intersect_{bed}_dm_fmin_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_{ref_for_ks_test}_oth_{next_bw}_eom.gplt"
    shell: 
        "sh scripts/plot_fsize.sh {input.enrichment_signal} {input.reassigned_bed_tsv}" 
        " {input.base_gnuplt_file} {output.eom_plt} {output.eom_eps} {output.eom_pdf}"
        " {wildcards.sample} {wildcards.slop} {wildcards.fsize} {input.colmean_eom}" 

rule savgol_smooth_other_bw: 
    input:
        #enrichment_signal = "other_bw_map/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_another_bw_{next_bw}.csv.gz_e.csv.gz",
        enrichment_signal = "other_bw_map/{sample}_intersect_{bed}_dm_fmin_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_{ref_for_ks_test}_oth_{next_bw}.csv.gz_e.csv.gz",
        
    params:
    output:
        savgol_smooth = "savgol_smooth_other_bw_manual/{sample}_intersect_{bed}_dm_fmin_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_{ref_for_ks_test}_oth_{next_bw}_sav_{window}_or_{order}.csv.gz_e.csv.gz"
    shell:
        "python $NGS_SCRIPTS_DIR/savgol_smooth.py {input.enrichment_signal} {output.savgol_smooth} {wildcards.window} {wildcards.order}"

def get_total_slop (wildcards):
    total_slop = 0 
    # get exact length of tfbs + slop 
    fp = open ("input_bed/"+ wildcards.bed + ".bed", "r") 
    l = fp.readline() 
    d_loc = [m.start() for m in re.finditer("\t", l)]   
    width = int(l[d_loc[1] + 1:d_loc[2]]) - int(l[d_loc[0] + 1:d_loc[1]]) 
    total_slop = width + int(wildcards.slop)
    return total_slop
    
rule colmean_savgol_other_bw:
    input:
        #savgol_smooth = "savgol_smooth_other_bw_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_another_bw_{next_bw}_savgol_{window}_order_{order}.csv.gz_e.csv.gz"
        savgol_smooth = "savgol_smooth_other_bw_manual/{sample}_intersect_{bed}_dm_fmin_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_{ref_for_ks_test}_oth_{next_bw}_sav_{window}_or_{order}.csv.gz_e.csv.gz"
    params:
        #total_slop = lambda wildcards: str(int(wildcards.slop) + int(wildcards.bed.split("bp_")[0]))
        total_slop = lambda wildcards: get_total_slop(wildcards)
    output:
        colmean_eom = "other_bw_colmean_savgol/{sample}_intersect_{bed}_dm_fmin_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_{ref_for_ks_test}_oth_{next_bw}_savgol_{window}_or_{order}_cm_eom.tsv" 
    shell:
        "sh scripts/calculate_colmeans_of_mnase_eom.sh {input.savgol_smooth}"
        " {wildcards.nclust} {params.total_slop} {output.colmean_eom}" 

rule colmean_plot_savgol_other_bw:
    input:
        #colmean_eom = "other_bw_colmean_savgol/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_another_bw_{next_bw}_savgol_{window}_order_{order}_colmean_eom.tsv", 
        colmean_eom = "other_bw_colmean_savgol/{sample}_intersect_{bed}_dm_fmin_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_{ref_for_ks_test}_oth_{next_bw}_savgol_{window}_or_{order}_cm_eom.tsv",
        colmean_line_colors = "utils/gnuplot_base_files/boxplot_colors.yaml" 
        
    params:
    output:
        colmean_eom_pdf = "plots/other_bw_colmean_savgol/{sample}_intersect_{bed}_dm_fmin_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_{ref_for_ks_test}_oth_{next_bw}_sav_{window}_or_{order}_cm_eom.pdf", 
        colmean_eom_png = "plots/other_bw_colmean_savgol/{sample}_intersect_{bed}_dm_fmin_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_{ref_for_ks_test}_oth_{next_bw}_sav_{window}_or_{order}_cm_eom.png", 

    shell:
        "Rscript scripts/heatmap_colmean_plot.R {input.colmean_eom}" 
        " {output.colmean_eom_pdf} {output.colmean_eom_png} {input.colmean_line_colors} {wildcards.nclust}"
rule other_bw_savgol_heatmaps:
    input:
        #savgol_smooth = "savgol_smooth_other_bw_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_another_bw_{next_bw}_savgol_{window}_order_{order}.csv.gz_e.csv.gz",
        savgol_smooth = "savgol_smooth_other_bw_manual/{sample}_intersect_{bed}_dm_fmin_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_{ref_for_ks_test}_oth_{next_bw}_sav_{window}_or_{order}.csv.gz_e.csv.gz",
        reassigned_bed_tsv = "fragment_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_reassigned_bed.tsv", # taken from heatmap_of_fragments_in_clsuters.smk
        base_gnuplt_file = lambda wildcards: get_gnuplot_base_file(wildcards),
        #colmean_eom = "other_bw_colmean_savgol/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_another_bw_{next_bw}_savgol_{window}_order_{order}_colmean_eom.tsv", 
        colmean_eom = "other_bw_colmean_savgol/{sample}_intersect_{bed}_dm_fmin_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_{ref_for_ks_test}_oth_{next_bw}_savgol_{window}_or_{order}_cm_eom.tsv" 
    params:
    output:
        eom_eps = "plots/other_bw_savgol_heatmaps/{sample}_intersect_{bed}_dm_fmin_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_{ref_for_ks_test}_oth_{next_bw}_sav_{window}_or_{order}.eps", 
        eom_pdf = "plots/other_bw_savgol_heatmaps/{sample}_intersect_{bed}_dm_fmin_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_{ref_for_ks_test}_oth_{next_bw}_sav_{window}_or_{order}.pdf", 
        eom_plt = "plots/other_bw_savgol_heatmaps/{sample}_intersect_{bed}_dm_fmin_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_{ref_for_ks_test}_oth_{next_bw}_sav_{window}_or_{order}.gplt", 
   
    shell: 
        "sh scripts/plot_fsize.sh {input.savgol_smooth} {input.reassigned_bed_tsv}" 
        " {input.base_gnuplt_file} {output.eom_plt} {output.eom_eps} {output.eom_pdf}"
        " {wildcards.sample} {wildcards.slop} {wildcards.fsize} {input.colmean_eom}"  
        

























############## Archived ################### 


#############################################################
##### Manually assign cluster id by recording appropriate ###
##### cluster id in file with visual inspection #############
##### its almost going to be a copy of boxplot rule above ###
#############################################################
#
#def get_manual_cl_id_file_bxplt(wildcards): 
#    # manually assinged cluster id can be found in `manual_cluster_assignment_files/`
#    # with file names as <sample>_<bed>_<min_flen>_<max_flen>_<min_reads>_<bw>_<nclust>_<max_iter>_cl_id.tsv
#    fname = "_".join([wildcards.sample, wildcards.bed, wildcards.min_flen,
#                      wildcards.max_flen, wildcards.min_reads, wildcards.bw,
#                      wildcards.nclust, wildcards.max_iter])
#    fname = "manual_cluster_assignment_files/" + fname + ".tsv"
#    return fname
#
#rule manual_box_plots_of_clusters:
#    input:
#        mean_line_plot_of_clusters = "plots/mean_plots_of_clusters/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_line_plot.png",
#        heatmap_file_eps = "plots/kmeans_kde_heatmap/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_heatmap.eps", 
#        exact_signal = "bw_map/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_chip_{chip}.csv.gz_exact_stats.tsv",
#        gnuplot_color_file = "utils/gnuplot_base_files/boxplot_colors.yaml",
#        #manual_cl_id_file = lambda wildcards: get_manual_cl_id_file (wildcards)
#        manual_cl_id_file = "manual_cluster_assignment_files/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_chip_{chip}_ref_cluster_{ref_for_ks_test}_manual_order.tsv",  #first get the order by chip and then use if for the manual plot  
#    params:
#        #system = config["system"], 
#        system = lambda wildcards: config["bed_annotation"][wildcards.bed[wildcards.bed.find("_")+1:len(wildcards.bed)]],
#        pseudo_chip = 0.01,
#	ref_for_ks_test = lambda wildcards: wildcards.ref_for_ks_test,
#        bxplt_min_value = config["plot_chip_params"]["min_value"],
#        bxplt_max_value = config["plot_chip_params"]["max_value"],
#        dotted_line_color = "\"#636363\"",
#        ref_ks_box_color = lambda wildcards: config["cluster_colors"][wildcards.ref_for_ks_test], 
#        source_annotation = lambda wildcards: config["source_annotation"][wildcards.chip], 
#        sample_title = lambda wildcards: config["sample_title"][wildcards.sample],
#        plot_title = lambda wildcards: get_boxplot_title(wildcards)
#    output:
#        chip_box_plot = "plots/manually_assigned_chip_boxplots/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_chip_{chip}_ref_cluster_{ref_for_ks_test}_bxplt_manual.png", 
#        chip_box_plot_pdf = "plots/manually_assigned_chip_boxplots/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_chip_{chip}_ref_cluster_{ref_for_ks_test}_bxplt_manual.pdf"
#    shell:
#        "sh scripts/plot_chip_manual.sh {input.exact_signal} {output.chip_box_plot}"
#        " {params.system} {params.pseudo_chip} {params.sample_title}"
#        " {wildcards.chip} \"{params.ref_for_ks_test}\" {params.bxplt_min_value}"
#        " {params.bxplt_max_value} {params.dotted_line_color} {params.ref_ks_box_color}"
#        " {input.gnuplot_color_file} \"{params.source_annotation}\""
#        " {output.chip_box_plot_pdf} \"{params.plot_title}\" {input.manual_cl_id_file}"  
#        
#
#
#
#def get_boxplot_title_sites_selected(wildcards):
#    tfbs_flank = "".join(["TFBS", u"\u00B1", wildcards.bed.split("_")[0]])
#    bw_str = "=".join(["bw", wildcards.bw])
#    title = "-".join([config["system"],
#                      config["sample_title"][wildcards.sample],
#                      bw_str, tfbs_flank, wildcards.site_string])
#    return title 
#rule plot_selected_chip_scores: 
#    input:
#        mean_line_plot_of_clusters = "plots/mean_plots_of_clusters/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_line_plot.png",
#        heatmap_file_eps = "plots/kmeans_kde_heatmap/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_heatmap.eps", 
#        exact_signal = "bw_map/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_chip_{chip}.csv.gz_exact_stats.tsv",
#        gnuplot_color_file = "utils/gnuplot_base_files/boxplot_colors.yaml"
#    params:
#        system = config["system"], 
#        pseudo_chip = 0.01,
#	ref_for_ks_test = lambda wildcards: wildcards.ref_for_ks_test,
#        bxplt_min_value = config["plot_chip_params"]["min_value"],
#        bxplt_max_value = config["plot_chip_params"]["max_value"],
#        dotted_line_color = "\"#636363\"",
#        ref_ks_box_color = lambda wildcards: config["cluster_colors"][wildcards.ref_for_ks_test], 
#        source_annotation = lambda wildcards: config["source_annotation"][wildcards.chip], 
#        sample_title = lambda wildcards: config["sample_title"][wildcards.sample],
#        plot_title = lambda wildcards: get_boxplot_title_sites_selected(wildcards)
#    output:
#        chip_box_plot = "plots/sites_selected_chip_box_plots/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_chip_{chip}_ref_cluster_{ref_for_ks_test}_bxplt_sites_{site_string}.png", 
#        chip_box_plot_pdf = "plots/sites_selected_chip_box_plots/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_chip_{chip}_ref_cluster_{ref_for_ks_test}_bxplt_sites_{site_string}.pdf",
#        ratio_barplot = "plots/sites_selected_chip_box_plots/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_chip_{chip}_ref_cluster_{ref_for_ks_test}_bxplt_sites_{site_string}.barplot.png"
#    shell:
#        "sh scripts/plot_chip_sites_selected.sh"
#        " {input.exact_signal} {output.chip_box_plot} {params.system}"
#        " {params.pseudo_chip} {params.sample_title} {wildcards.chip}"
#        " \"{params.ref_for_ks_test}\" {params.bxplt_min_value}"
#        " {params.bxplt_max_value} {params.dotted_line_color}"
#        " {params.ref_ks_box_color} {input.gnuplot_color_file}"
#        " \"{params.source_annotation}\" {output.chip_box_plot_pdf}"
#        " \"{params.plot_title}\" \"{wildcards.site_string}\""
#        " {output.ratio_barplot}"








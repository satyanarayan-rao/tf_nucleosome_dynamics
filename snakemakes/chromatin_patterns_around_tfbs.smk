rule get_mnase_scores: 
    input:
        bed_file = "input_bed/{bed}.bed", 
        bw_file = lambda wildcards: config["chip_bigwigs"][wildcards.chip]
    params:
    output:
        mnase_score = "mnase_score/{chip}_on_{bed}.csv.gz",
        mnase_eom = "mnase_score/{chip}_on_{bed}.csv.gz_e.csv.gz", 
        mnase_exact_stats = "mnase_score/{chip}_on_{bed}.csv.gz_exact_stats.tsv" 
    shell:
        "python $NGS_SCRIPTS_DIR/map_bw_to_bed.py {input.bw_file}"
        " {input.bed_file} {output.mnase_score}"
rule eom_to_eps: 
    input:
        mnase_eom = "mnase_score/{chip}_on_{bed}.csv.gz_e.csv.gz",
        gplt_file = lambda wildcards: config["gnuplt_prefix"][wildcards.size_label], 
        bed_file = "input_bed/{bed}.bed"
    params: 
        xlim = lambda wildcards: wildcards.bed.split ("bp_")[0]
    output:
        mnase_eps = "plots/eom_eps_pdf/heatmap_{chip}_on_{bed}_eom_{size_label}.eps", 
        mnase_pdf = "plots/eom_eps_pdf/heatmap_{chip}_on_{bed}_eom_{size_label}.pdf", 
        mnase_gplt = "plots/eom_eps_pdf/heatmap_{chip}_on_{bed}_eom_{size_label}.gplt"
        
    shell:
        "sh scripts/matrix_gz_to_plt.sh {input.mnase_eom} {input.gplt_file}"
        " {output.mnase_eps} {output.mnase_pdf} {output.mnase_gplt}"
        " {wildcards.chip} {input.bed_file} {params.xlim}" 
rule mean_eom_plots: 
    input:
        mnase_eom = "mnase_score/{chip}_on_{bed}.csv.gz_e.csv.gz",
    params:
    output:
        colmeans = "mnase_score/{chip}_on_{bed}_eom_colmeans.tsv"
    shell:
        "python $NGS_SCRIPTS_DIR/colmeans_gz_file_with_row_names.py"
        " {input.mnase_eom} {output.colmeans}"

rule plot_colmeans: 
    input: 
        files = ["mnase_score/mnase_mcf7_on_750bp_chip_positive_eom_colmeans.tsv",
                 "mnase_score/mnase_mcf7_on_750bp_chip_negative_eom_colmeans.tsv"]
    params:
        
    output:
        cmean_plot_eps = "plots/eom_eps_pdf/positive_and_negative.eps", 
        cmean_plot_pdf = "plots/eom_eps_pdf/positive_and_negative.pdf", 
        cmean_plot_gplt = "plots/eom_eps_pdf/positive_and_negative.gplt" 
    shell:
        "sh scripts/colmeans_plt.sh \"{input.files}\" {output.cmean_plot_eps}"
        " {output.cmean_plot_pdf} {output.cmean_plot_gplt}"



rule kmeans_on_mnase: 
    input:
        mnase_eom = "mnase_score/{chip}_on_{bed}.csv.gz_e.csv.gz", 
    params:
        
    output:
        kmeans_out = "kmeans_on_mnase_scores/kmeans_{chip}_on_{bed}_nclust_{nclust}.csv.gz", 
        row_order = "kmeans_on_mnase_scores/kmeans_{chip}_on_{bed}_nclust_{nclust}_row_order.tsv", 
        kmeans_rdata = "kmeans_on_mnase_scores/kmeans_{chip}_on_{bed}_nclust_{nclust}.Rdata",
    shell:
        "Rscript scripts/kmeans_cluster_eq_n.R {input.mnase_eom}"
        " {output.kmeans_out} {output.row_order} {wildcards.nclust}"
        " {output.kmeans_rdata}"
def get_kmeans_rdata_files (wildcards):
    flist = []
    for cl in config["elbow_plot_kmeans"][wildcards.bed]["cl_set"]:
        flist.append("kmeans_on_mnase_scores/kmeans_{chip}_on_{bed}_nclust_{nclust}.Rdata".format(
                      chip = wildcards.chip, nclust = cl, bed = wildcards.bed))
    
    return flist

rule elbow_plot_kmeans: 
    input: 
        rdata_files = lambda wildcards: get_kmeans_rdata_files(wildcards) 
    params:  
        cl_ids = lambda wildcards: "@".join(map(str,
                     config["elbow_plot_kmeans"][wildcards.bed]["cl_set"]))
    output:
        elbow_tsv = "kmeans_on_mnase_scores/kmeans_{chip}_on_{bed}.elbow.tsv", 
        elbow_png = "plots/kmeans_on_mnase_eom/kmeans_{chip}_on_{bed}.elbow.png", 
        
    shell:
        "Rscript scripts/elbow_chip_kmeans.R \"{input.rdata_files}\" {params.cl_ids} {output.elbow_tsv}"
        " {output.elbow_png} {wildcards.chip}"
 
        
rule plot_kmeans_mnase_eom:
    input: 
        kmeans_out = "kmeans_on_mnase_scores/kmeans_{chip}_on_{bed}_nclust_{nclust}.csv.gz", 
        row_order = "kmeans_on_mnase_scores/kmeans_{chip}_on_{bed}_nclust_{nclust}_row_order.tsv", 
        #gnuplt_base_file = "utils/gnuplot_base_files/kmeans_dist.gplt"
        gnuplt_base_file = "utils/gnuplot_base_files/fragment_matrix_130_180.gplt"
    params:
        label = "Kmeans-Nuc", 
        slop = lambda wildcards: wildcards.bed.split("bp_")[0]
    output:
        kmeans_gnuplt_file = "plots/kmeans_on_mnase_eom/kmeans_{chip}_on_{bed}_nclust_{nclust}.gplt", 
        kmeans_plot_eps = "plots/kmeans_on_mnase_eom/kmeans_{chip}_on_{bed}_nclust_{nclust}.eps",
        kmeans_plot_pdf = "plots/kmeans_on_mnase_eom/kmeans_{chip}_on_{bed}_nclust_{nclust}.pdf"
    shell:
        "sh scripts/plot_kmeans.sh {input.kmeans_out} {input.gnuplt_base_file}"
        " {output.kmeans_gnuplt_file} {output.kmeans_plot_eps} {output.kmeans_plot_pdf}"
        " {params.slop} {params.slop} {input.row_order} {params.slop} \"Nuc\" \"Nuc\" \"{params.label}\""




rule do_ccr_mnase: 
    input:
        mnase_score = "mnase_score/{chip}_on_{bed}.csv.gz"
    params:
    output:
        ccr_mnase = "mnase_ccr/ccr_{chip}_on_{bed}_lag_{lag}.csv.gz", 
        nan_mnase = "mnase_ccr/nan_{chip}_on_{bed}_lag_{lag}.csv.gz"
    shell:
        "Rscript scripts/cross_correlation_a_and_b.R"
        " {input.mnase_score} {input.mnase_score}"
        " {wildcards.lag} {output.ccr_mnase} {output.nan_mnase}"   

rule do_kmeans_on_ccr_mnase:
    input:
        ccr_mnase = "mnase_ccr/ccr_{chip}_on_{bed}_lag_{lag}.csv.gz"
    params:
    output:
        kmeans_mnase = "kmeans_mnase/kmeans_{chip}_on_{bed}_lag_{lag}_nclust_{nclust}.csv.gz", 
        row_order = "kmeans_mnase/kmeans_{chip}_on_{bed}_lag_{lag}_nclust_{nclust}_row_order.tsv"
    shell:
        "Rscript scripts/kmeans_cluster_eq_n.R {input.ccr_mnase}"
        " {output.kmeans_mnase} {output.row_order} {wildcards.nclust}"
rule plot_kmeans_mnase: 
    input:
        kmeans_mnase = "kmeans_mnase/kmeans_{chip}_on_{bed}_lag_{lag}_nclust_{nclust}.csv.gz", 
        row_order = "kmeans_mnase/kmeans_{chip}_on_{bed}_lag_{lag}_nclust_{nclust}_row_order.tsv", 
        gnuplt_base_file = "utils/gnuplot_base_files/kmeans_dist.gplt"
    params:
        ccr_label = "AutoCorrelation", 
        slop = lambda wildcards: wildcards.bed.split("bp_")[0]
    output:
        kmeans_gnuplt_file = "plots/kmeans_on_mnase/kmeans_{chip}_on_{bed}_lag_{lag}_nclust_{nclust}.gplt", 
        kmeans_plot_eps = "plots/kmeans_on_mnase/kmeans_{chip}_on_{bed}_lag_{lag}_nclust_{nclust}.eps",
        kmeans_plot_pdf = "plots/kmeans_on_mnase/kmeans_{chip}_on_{bed}_lag_{lag}_nclust_{nclust}.pdf",
    shell:
        "sh scripts/plot_kmeans.sh {input.kmeans_mnase} {input.gnuplt_base_file} {output.kmeans_gnuplt_file} {output.kmeans_plot_eps} {output.kmeans_plot_pdf} {params.slop} {params.slop} {input.row_order} {params.slop} \"Nuc\" \"Nuc\" \"{params.ccr_label}\""
        
      

#rule get_mnase_scores: 
#    input:
#    params:
#    output:
#    shell:


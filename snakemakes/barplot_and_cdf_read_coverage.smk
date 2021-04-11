rule prepare_read_coverage_per_site: 
    input:
        cnt_matrix = "flen_count_matrix/{sample}_intersect_{bed}_flen_in_motif.tsv", 
        starting_bed = "input_bed/{bed}.bed"
    params:
    output:
        barplot_tsv = "flen_barplot_data/{sample}_intersect_{bed}_barplot_data.tsv"
    shell:
        "python  scripts/prepare_barplot_data.py {input.cnt_matrix}"
        " {input.starting_bed} {output.barplot_tsv}" 
rule plot_postive_vs_negative_hist: 
    input:
        barplot_tsv = "flen_barplot_data/{sample}_intersect_{bed}_barplot_data.tsv"
    params:
    output:
        hist_plot = "plots/cnt_hist/{sample}_intersect_{bed}_barplot_data.hist.png", 
        cdf_plot = "plots/cnt_hist/{sample}_intersect_{bed}_barplot_data.cdf.png", 
        hist_tsv = "flen_barplot_data/{sample}_intersect_{bed}_barplot_data.hist.tsv",
        cdf_tsv = "flen_barplot_data/{sample}_intersect_{bed}_barplot_data.cdf.tsv"
        
    shell:
  #      "Rscript scripts/cnt_hist.R {input.barplot_tsv} {output.hist_plot}"
        "Rscript scripts/cnt_hist_capped.R {input.barplot_tsv} {output.hist_plot}"
        " {output.cdf_plot} {wildcards.sample} {output.hist_tsv} {output.cdf_tsv}"  

rule do_kde_of_flen: 
    input:
        flen_count_in_moitf = "flen_count_matrix/{sample}_intersect_{bed}_flen_in_motif.tsv", 
        genomewide_background = lambda wildcards: config["genomic_len_dist"][wildcards.sample]
    params:
        density_n = 100 # number of points where KDE should be calcuated
    output:
        flen_kde_matrix = "kde_matrix/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}.tsv"
    shell:
        "Rscript scripts/density_of_flen_in_motif.R {input.flen_count_in_moitf} {params.density_n} {wildcards.min_flen} {wildcards.max_flen} {output.flen_kde_matrix} {wildcards.min_reads} {wildcards.bw} {input.genomewide_background}"
    
rule kmeans_on_kde:
    input:
        flen_kde_matrix = "kde_matrix/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}.tsv"
    params:
    output:
        out_kmeans_file = "kmeans_on_kde/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}.tsv", 
        out_kmeans_file_bed_tsv = "kmeans_on_kde/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_bed.tsv"
        
    shell:
        "sh scripts/kmeans_on_kde.sh {input.flen_kde_matrix} {wildcards.nclust} {wildcards.max_iter} {output.out_kmeans_file} {output.out_kmeans_file_bed_tsv}"

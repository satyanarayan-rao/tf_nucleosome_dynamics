rule count_fsize_from_bedgz: 
    input:
        input_file = lambda wildcards: config["bedgz_annotation"][wildcards.sample]
    params:
    output:
        count_file = "overall_flen_distribution/flen_dist_{sample}.count",
        hist_file = "overall_flen_distribution/flen_dist_{sample}.hist"
    shell:
        "python scripts/bedgz2flen_hist.py {input.input_file}"
        " {output.hist_file} {output.count_file}" 

def get_hist_files(wildcards):
    flist = []
    for f in config["overall_flen_dist"][wildcards.setting]["systems"]:
        flist.append("overall_flen_distribution/flen_dist_{s}.hist".format(s=f)) 
    return flist
def get_count_files(wildcards):
    flist = []
    for f in config["overall_flen_dist"][wildcards.setting]["systems"]:
        flist.append("overall_flen_distribution/flen_dist_{s}.count".format(s=f)) 
    return flist
rule plot_len_dist_multiple_samples: 
    input:
        input_hist_files = lambda wildcards: get_hist_files(wildcards),
        input_cnt_files = lambda wildcards: get_count_files(wildcards)
    params:
        labels = lambda wildcards: config["overall_flen_dist"][wildcards.setting]["labels"]
    output:
        len_dist_multiple_sample = "plots/overall_len_dist/setting_{setting}_plot.png"
    shell:
        "Rscript scripts/multiple_sample_hist_plot.R \"{input.input_hist_files}\""
        " \"{input.input_cnt_files}\" \"{params.labels}\" {output.len_dist_multiple_sample}"

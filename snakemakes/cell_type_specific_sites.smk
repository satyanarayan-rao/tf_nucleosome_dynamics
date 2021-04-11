rule get_chip_scores: 
    input:
        input_bed = "input_bed/{bed}.bed", 
        bw_file = lambda wildcards: config["chip_bigwigs"][wildcards.chip]
    params:
    output: 
        chip_scores = "cell_type_specific_sites/{chip}_on_{bed}_exact_stats.tsv" 
    shell:
        "python $NGS_SCRIPTS_DIR/map_bw_exact_stats_to_bed.py {input.bw_file}"
        " {input.input_bed} {output.chip_scores}"

def get_files_for_kmeans (wildcards):
    file_list = config["chip_kmeans"][wildcards.setting]["file_list"]
rule kmeans_on_chip_scores: 
    input:
        inp_files = lambda wildcards: config["chip_kmeans"][wildcards.setting]["file_list"]
    params:
        label_list = lambda wildcards: "@".join(config["chip_kmeans"][wildcards.setting]["labels"])
    output: 
        kmeans_out_file = "chip_kmeans/setting_{setting}_chip_kmeans.tsv" 
    shell:
        "Rscript scripts/chip_kmeans.R \"{input.inp_files}\" \"{params.label_list}\"" 
        " {output.kmeans_out_file}" 








#rule get_chip_scores: 
#    input:
#    params:
#    output: 
#    shell:
#     

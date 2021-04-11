rule extract_disease_subtype_count_matrix: 
    input:
        raw_count_file = "/beevol/home/satyanarr/data/data_from_papers/corces_mr_et_al_science_2018/cancer_type_specific_count_matrices/BRCA_raw_counts.txt", 
        patient_id_file = lambda wildcards: config["tcga_atac_params"][wildcards.disease_subtype]
    params:
    output:
        disease_subtype_count_matrix = "tcga_atac_seq_analysis/raw_count_matrix_{disease_subtype}.tsv"
    shell:
        "sh scripts/exctact_disease_subtype_columns.sh {input.raw_count_file} {input.patient_id_file} {output.disease_subtype_count_matrix}"

rule intesect_binding_sites_with_count_matrix: 
    input:
        disease_subtype_count_matrix = "tcga_atac_seq_analysis/raw_count_matrix_{disease_subtype}.tsv", 
        bed_file = lambda wildcards: config["bedfile_annotation"][wildcards.bed]
    params: 
        chunk_size = 10000
    output: 
        tfbs_mapped_to_count_matrix = "tcga_atac_seq_analysis/tfbs_{bed}_mapped_to_{disease_subtype}.tsv"
    shell:
        "sh scripts/intersect_tfbs_with_count_matrix.sh {input.bed_file}"
        " {input.disease_subtype_count_matrix} {output.tfbs_mapped_to_count_matrix}" 
rule prepare_file_for_boxplot:
    input: 
        tfbs_mapped_to_count_matrix = "tcga_atac_seq_analysis/tfbs_{bed}_mapped_to_{disease_subtype}.tsv"
    params:
    output:
        long_listed_tsv_with_class_label = "tcga_atac_seq_analysis/boxplot_data_for_{bed}_mapped_to_{disease_subtype}.tsv"
    shell:
        "python scripts/preapre_boxplot_data.py {input.tfbs_mapped_to_count_matrix}"
        " {output.long_listed_tsv_with_class_label} "

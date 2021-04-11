rule append_cell_line_in_bed: 
    input:
        motif_file = "/beevol/home/satyanarr/workplace/projects/cfDNA-TF-footprints/TFs/FOXA1/ChIP/WgEncodeRegTfbsClustered/V3/cfDNA_to_ChIP_motif_distance/fimo_on_foxa1/fimo.bed", 
        cell_line_pkl = "/beevol/home/satyanarr/workplace/projects/cfDNA-TF-footprints/TFs/FOXA1/ChIP/WgEncodeRegTfbsClustered/V3/cfDNA_to_ChIP_motif_distance/fimo_on_foxa1/fimo.bedm2c.pkl"
    params:
    output:
        motif_with_cell_line = "input_bed/encode_foxa1_motifs_with_cell_line.bed"
    shell:
        "python scripts/encode_append_cell_line.py {input.motif_file} {input.cell_line_pkl}"
        " {output.motif_with_cell_line}"

def get_tissue_cell_lines():
    tissue_cell_lines = ""
    tissue_list = []
    for k in config["tissue_type_based_sites"]["encode"]:
        c_string = config["tissue_type_based_sites"]["encode"][k]
        tissue_list.append("@".join([k, c_string])) 
    tissue_cell_lines = "^".join(tissue_list)
    return tissue_cell_lines

rule tissue_based_bed_mututally_exclusive: 
    input: 
        motif_with_cell_line = "input_bed/encode_foxa1_motifs_with_cell_line.bed"
    params:
        tissue_cell_lines = get_tissue_cell_lines()
    output:
        tissue_bed = "input_bed/encode_tissue_based_motifs.bed" 
    shell:
        "python scripts/encode_append_tissue.py {input.motif_with_cell_line} \"{params.tissue_cell_lines}\""
        " {output.tissue_bed}"

rule append_tissue_in_carroll_sites: 
    input: 
        carroll_fimo_file = "/beevol/home/satyanarr/workplace/projects/cfDNA-TF-footprints/TFs/FOXA1/ChIP/Carroll_JS_Lab/new_cfDNA_to_ChIP_motif_distance/fimo_on_foxa1/fimo.bed"
    params:
    output:
        appended_tissue_type = "input_bed/carroll_tissue_based_motifs.bed"
    shell:
        "awk '{{print $1\"\t\"$2\"\t\"$3\"\t\"$4\"`BreastChIP_Carroll\t\"$5\"\t\"$6}}'"
        " {input.carroll_fimo_file} >  {output.appended_tissue_type}"  
rule append_neg_ctl_label_in_segway: 
    input:
        neg_ctl_segway_bed = "input_bed/neg_ctl_segway_jaspar_MA0148.3_num_sites_selected_577065.bed"
    params:
    output:
        neg_ctl_segway_with_label = "input_bed/neg_ctl_segway_jaspar_MA0148.3_num_sites_selected_577065_with_label.bed"
    shell:
        "awk '{{print $1\"\t\"$2\"\t\"$3\"\t\"$4\"`neg_ctl\t\"$5\"\t\"$6}}'"
        " {input.neg_ctl_segway_bed} >  {output.neg_ctl_segway_with_label}"  
    

rule annotate_encode_carroll_commons:
    input:
        encode_file = "input_bed/encode_tissue_based_motifs.bed", # look at the rule tissue_bed `tissue_based_bed_mututally_exclusive` = "input_bed/encode_tissue_based_motifs.bed" 
        carroll_file = "input_bed/carroll_tissue_based_motifs.bed", # look at the rule `append_tissue_in_carroll_sites`;  appended_tissue_type = "input_bed/carroll_tissue_based_motifs.bed" 
        accepted_chrom = "metadata/accepted_chrom.tsv"
    params:
    output:
        unique_or_common_sites = "input_bed/unique_encode_and_carroll_motifs.bed", 
        accepted_unique_or_common_sites = "input_bed/unique_encode_and_carroll_motifs_accepted_chrom.bed"
    shell:
        "python scripts/assign_encode_or_carroll_label.py {input.encode_file}"
        " {input.carroll_file} {output.unique_or_common_sites}"
        "; python $NGS_SCRIPTS_DIR/select_for_accepted_chr.py"
        " --start_bed {output.unique_or_common_sites}"
        " --accepted_chrom {input.accepted_chrom}"
        " --out_bed {output.accepted_unique_or_common_sites}" 

rule append_encode_carroll_and_neg_ctl: 
    input:
        accepted_unique_or_common_sites = "input_bed/unique_encode_and_carroll_motifs_accepted_chrom.bed", 
        neg_ctl_segway_with_label = "input_bed/neg_ctl_segway_jaspar_MA0148.3_num_sites_selected_577065_with_label.bed"
    params:
    output:
        appended_bed = "input_bed/encode_carroll_and_neg_ctl.bed"
    shell: 
        "cat {input.accepted_unique_or_common_sites} {input.neg_ctl_segway_with_label} > {output.appended_bed}"
rule slop_accepted_bed: 
    input:
        appended_bed = "input_bed/encode_carroll_and_neg_ctl.bed",
        hg38_genome = "metadata/hg38.genome"
    params:
    output:
        slopped_bed = "input_bed/encode_carroll_and_neg_ctl_slop_{slop}.bed" 
    shell: 
        "bedtools slop -b {wildcards.slop} -g {input.hg38_genome} -i"
        " {input.appended_bed} > {output.slopped_bed}"
rule map_exact_stats: 
    input:
        slopped_bed = "input_bed/encode_carroll_and_neg_ctl_slop_{slop}.bed" ,
        bw_file = lambda wildcards: config["chip_bigwigs"][wildcards.chip_bigwig]
    params:
    output:
        exact_stats_file = "chip_analysis_for_site_selection/bw_{chip_bigwig}_mapped_to_subset_{tissue}_exact_stats_for_slop_{slop}.tsv"
    shell:
        "sh scripts/chip_exact_stats.sh {input.slopped_bed} \"{wildcards.tissue}\""
        " {input.bw_file}  {output.exact_stats_file}"

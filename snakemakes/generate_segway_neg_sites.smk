def prepare_annotation_string(wildcards):
    annot_vec = config["genome_annotations"]["segway"][wildcards.regions]
    annot_str = "|".join(annot_vec)
    return annot_str
rule liftover_hg19_to_hg38: 
    input: 
        bed_file = "input_bed/MCF-7.bed", 
        chain_file = "/beevol/home/satyanarr/data/ucsc/liftover_chain_files/hg19ToHg38.over.chain"
        
    params: 
    output:
        hg38_bed_file = "input_bed/MCF-7-hg38.bed", 
        unmapped = "input_bed/MCF-7-hg38-unmapped.bed"
    shell:
        "liftOver {input.bed_file} {input.chain_file}" 
        " {output.hg38_bed_file} {output.unmapped}"

rule select_defined_region: 
    input: 
        hg38_bed_file = "input_bed/MCF-7-hg38.bed", 
    params:
        annotation_in_region = lambda wildcards: prepare_annotation_string(wildcards)
    output: 
        selected = "input_bed/selected_mcf7_{regions}.bed"
    shell:
        "egrep \"{params.annotation_in_region}\" {input.hg38_bed_file} > {output.selected}"

rule prepare_formatted_bed:
    input: 
        selected = "input_bed/selected_mcf7_{regions}.bed"
    params:
    output: 
        formatted = "input_bed/formatted_mcf7_{regions}.bed"
    shell: 
        "python $NGS_SCRIPTS_DIR/format_bed.py {input.selected} {output.formatted}"
    
rule get_fasta: 
    input: 
        formatted = "input_bed/formatted_mcf7_{regions}.bed", 
        hg38_fa = "/beevol/home/satyanarr/data/ucsc/human_genome/hg38/hg38.fa"
    params:
    output:
        formatted_fasta = "input_bed/fasta_formatted_mcf7_{regions}.fa"
    shell:
        "bedtools getfasta -fi {input.hg38_fa} -fo {output.formatted_fasta}"
        " -bed {input.formatted} -name" 

rule run_fimo: 
    input: 
        formatted_fasta = "input_bed/fasta_formatted_mcf7_{regions}.fa",
        motif_file = lambda wildcards: config["motif"][wildcards.motif]
    params:
        dest_dir = lambda wildcards: "fimo/fimo_on_segway_" + wildcards.regions
    output: 
        fimo_out_segway = "fimo/fimo_on_segway_{regions}/fimo_{motif}.txt"
    shell:
        "sh scripts/run_fimo.sh {input.formatted_fasta} {input.motif_file}"
        " {params.dest_dir} {output.fimo_out_segway}"

rule select_for_unmasked:  # select for discovered motifs with all capital letters
    input: 
        fimo_out_segway = "fimo/fimo_on_segway_{regions}/fimo_{motif}.txt"
    params:
    output: 
        unmasked_motifs = "fimo/fimo_on_segway_{regions}/unmasked_fimo_{motif}.txt"
    shell: 
        "head -1 {input.fimo_out_segway} > {output.unmasked_motifs} ; awk '{{if (toupper($NF) == $NF) {{print $0}} }}' {input.fimo_out_segway} >> {output.unmasked_motifs}" 

rule txt_to_bed: 
    input:
        unmasked_motifs = "fimo/fimo_on_segway_{regions}/unmasked_fimo_{motif}.txt"
    params:
    output: 
        bed = "fimo/fimo_on_segway_{regions}/unmasked_fimo_{motif}.bed"
    shell: 
        "python $NGS_SCRIPTS_DIR/fimo2bed.py --fimo {input.unmasked_motifs} --out {output.bed} "
rule bed_keep_accepted_chroms: 
    input: 
        bed = "fimo/fimo_on_segway_{regions}/unmasked_fimo_{motif}.bed", 
        accepted_chrom = "metadata/accepted_chrom.tsv"
    params: 
    output:
        accepted_bed = "fimo/fimo_on_segway_{regions}/accepted_chrom_unmasked_fimo_{motif}.bed" 
    shell: 
        "python $NGS_SCRIPTS_DIR/select_for_accepted_chr.py --start_bed {input.bed}"
        " --accepted_chrom {input.accepted_chrom} --out_bed {output.accepted_bed}"
rule rm_chip_motifs: 
    input:
        accepted_bed = "fimo/fimo_on_segway_{regions}/accepted_chrom_unmasked_fimo_{motif}.bed",
        chip_motif_pkl = "input_bed/unique_encode_and_carroll_with_strand_chrom_filtered_motif.pkl"
     
    params:
    output:
        chip_removed_bed = "fimo/fimo_on_segway_{regions}/chip_removed_fimo_{motif}.bed"
    shell:
        "python scripts/segway_remove_chip_sites_from_neg_ctl.py"
        " {input.accepted_bed} {input.chip_motif_pkl} {output.chip_removed_bed}"
rule segway_randomly_select_sites: 
    input: 
        chip_removed_bed = "fimo/fimo_on_segway_{regions}/chip_removed_fimo_{motif}.bed"
    params:
    output:
        shuf_line_file = "fimo/fimo_on_segway_{regions}/shuffled_line_id_{motif}_num_sites_selected_{sites}.tsv",
        shuffled_selected = "fimo/fimo_on_segway_{regions}/neg_ctl_segway_{motif}_num_sites_selected_{sites}.bed" 
    shell:
      "sh scripts/shuffle_using_seed.sh {input.chip_removed_bed} {output.shuf_line_file}"
      " {output.shuffled_selected} {wildcards.sites}" 

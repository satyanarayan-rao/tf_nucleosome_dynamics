rule foxa1_motif_in_hg38:
    input: 
        genome = "/beevol/home/satyanarr/data/ucsc/human_genome/hg38/hg38.fa", 
        motif_file = lambda wildcards: config["motif"][wildcards.motif]
    params:
    output:
        fimo_out_file = "fimo_for_neg_ctl/fimo_{motif}.txt" 
    shell:
        "sh scripts/run_fimo.sh {input.genome} {input.motif_file} fimo_for_neg_ctl {output.fimo_out_file}"
rule fimo_txt_to_fimo_bed:
    input:
        fimo_out_file = "fimo_for_neg_ctl/fimo_{motif}.txt" 
    params:
    output:
        fimo_bed = "fimo_for_neg_ctl/fimo_{motif}.bed"
    shell:
        "awk 'NR>1{{print $3\"\t\"$4\"\t\"$5\"\t\"$3\":\"$4\"-\"$5\"|\"$2\"`\"NR\"_\"$NF\"%neg_ctl\t\"$7\"\t\"$6}}' {input.fimo_out_file} > {output.fimo_bed}"

rule keep_accepted_chroms: 
    input:
        fimo_bed = "fimo_for_neg_ctl/fimo_{motif}.bed", 
        accepted_chrom = "metadata/accepted_chrom.tsv"
    params:
    output:
        accepted_bed = "fimo_for_neg_ctl/accepted_chrom_fimo_{motif}.bed"
    shell:
        "python $NGS_SCRIPTS_DIR/select_for_accepted_chr.py --start_bed {input.fimo_bed}"
        " --accepted_chrom {input.accepted_chrom}"
        " --out_bed {output.accepted_bed}"

rule create_motif_dict: 
    input:
        chip_bed = "input_bed/unique_encode_and_carroll_with_strand_chrom_filtered.bed"
    params:
    output:
        chip_motif_pkl = "input_bed/unique_encode_and_carroll_with_strand_chrom_filtered_motif.pkl"
    shell:
        "python scripts/build_motif_pkl.py {input.chip_bed} {output.chip_motif_pkl}"
rule remove_chip_motifs:
    input:
        accepted_bed = "fimo_for_neg_ctl/accepted_chrom_fimo_{motif}.bed", 
        chip_motif_pkl = "input_bed/unique_encode_and_carroll_with_strand_chrom_filtered_motif.pkl"
    params:
    output:
        chip_removed_bed = "fimo_for_neg_ctl/chip_removed_fimo_{motif}.bed"
    shell:
        "python scripts/remove_chip_sites_from_neg_ctl.py"
        " {input.accepted_bed} {input.chip_motif_pkl} {output.chip_removed_bed}"

rule randomly_select_sites: 
    input: 
        chip_removed_bed = "fimo_for_neg_ctl/chip_removed_fimo_{motif}.bed"
    params:
    output:
        shuf_line_file = "fimo_for_neg_ctl/shuffled_line_id_{motif}_num_sites_selected_{sites}.tsv",
        shuffled_selected = "fimo_for_neg_ctl/neg_ctl_fimo_{motif}_num_sites_selected_{sites}.bed"
    shell:
      "sh scripts/shuffle_using_seed.sh {input.chip_removed_bed} {output.shuf_line_file}"
      " {output.shuffled_selected} {wildcards.sites}"    

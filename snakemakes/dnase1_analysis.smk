
rule prepare_input_bed:
    input: 
        inp_bed = "input_bed/unique_encode_and_carroll_with_strand_chrom_filtered_no_chrY.bed"
    output:
        out_bed = "input_bed/dnase1_bed/chip_peaks_width_{width}.bed" 
    shell:
        "sh scripts/select_chip_peaks.sh {input.inp_bed} {wildcards.width} {output.out_bed}"
 
rule slop_input_bed: 
    input:
        out_bed = "input_bed/dnase1_bed/chip_peaks_width_{width}.bed", 
        genome_size = "metadata/hg38.genome"
    output:
        slop_bed = "input_bed/dnase1_bed/slop_{slop}_chip_peaks_width_{width}.bed" 
    shell:
        "bedtools slop -b {wildcards.slop} -g {input.genome_size} -i {input.out_bed} > {output.slop_bed}" 
rule map_dnase_to_bed:
    input:
        bigwig_file = "/beevol/home/satyanarr/workplace/data/geo/GSE51915_He_et_al_Nature_Methods_2013/{bigwig}.bigwig",
        slop_bed = "input_bed/dnase1_bed/slop_{slop}_chip_peaks_width_{width}.bed" 
    params:
    output:
        mapped_bedgz = "dnase1_analysis/{bigwig}_mapped_to_slop_{slop}_chip_peaks_width_{width}.csv.gz", 
        mapped_eom = "dnase1_analysis/{bigwig}_mapped_to_slop_{slop}_chip_peaks_width_{width}.csv.gz_e.csv.gz", 
        exact_signal = "dnase1_analysis/{bigwig}_mapped_to_slop_{slop}_chip_peaks_width_{width}.csv.gz_exact_stats.tsv"
    shell:
        "python $NGS_SCRIPTS_DIR/map_bw_to_bed.py {input.bigwig_file}"
        " {wildcards.slop} {input.bed_file} {output.mapped_bedgz}"

rule order_dnase_peaks: 
    input:
        mapped_eom = "dnase1_analysis/{bigwig}_mapped_to_slop_{slop}_chip_peaks_width_{width}.csv.gz_e.csv.gz", 
    params: 
    output:
        ordered_sites = "dnase1_analysis/dnase1_ordered_{bigwig}_mapped_to_slop_{slop}_chip_peaks_width_{width}_e.csv.gz", 
        bed_with_dnase_peak = "dnase1_analysis/nearest_dnase1_peak_{bigwig}_mapped_to_slop_{slop}_chip_peaks_width_{width}.bed"
    shell:
        "python scripts/call_dnase1_peak.py {input.mapped_eom}"
        " {wildcards.slop} {output.ordered_sites} {output.bed_with_dnase_peak}"
rule order_by_in_and_out: 
    input:
        bed_with_dnase_peak = "dnase1_analysis/nearest_dnase1_peak_{bigwig}_mapped_to_slop_{slop}_chip_peaks_width_{width}.bed"
        
    output:
        ordered_bed = "dnase1_analysis/in_and_out_{bigwig}_mapped_to_slop_{slop}_chip_peaks_width_{width}.bed" 
    shell:
        "sh scripts/order_by_in_and_out.sh {input.bed_with_dnase_peak}"
        " {output.ordered_bed}" 
#rule dnase_enrichment_matrix: 
#     input:
#        ordered_bed = "dnase1_analysis/in_and_out_{bigwig}_mapped_to_slop_{slop}_chip_peaks_width_{width}.bed", 
#        bigwig_file = "/beevol/home/satyanarr/workplace/data/geo/GSE51915_He_et_al_Nature_Methods_2013/{bigwig}.bigwig",

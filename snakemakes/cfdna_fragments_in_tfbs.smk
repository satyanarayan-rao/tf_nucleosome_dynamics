############### These two commented rules are no longer required as the third rule serves the purpose ############# 
#rule foxa1_bed_intersect_with_bedgz: 
#    input:
#        tfbs = "input_bed/{bed}.bed",
#        bedgz_file = lambda wildcards: config["bedgz_annotation"][wildcards.sample]
#    params:
#    output:
#        #bed_intersect_file = "bed_intersect/{sample}_intersect_{bed}.bed.gz"
#        bed_intersect_file = "bed_intersect/{sample}_intersect_{bed}.bed.gz"
#    shell:
#        "sh scripts/bed_intersect.sh {input.bedgz_file} {input.tfbs} {output.bed_intersect_file}" 
#rule count_flen_in_tfbs: 
#    input:
#        bed_intersect_file = "bed_intersect/{sample}_intersect_{bed}.bed.gz"
#    params:
#    output:
#        count_matrix = "flen_count_matrix/{sample}_intersect_{bed}_flen_count.tsv", 
#        flen_count_in_moitf = "flen_count_matrix/{sample}_intersect_{bed}_flen_in_motif.tsv"
#    shell:
#        "python scripts/generate_flen_count_matrix.py {input.bed_intersect_file} {output.count_matrix} {output.flen_count_in_moitf}"
#rule count_flen_in_tfbs_new:
#    input: 
#        tfbs = "input_bed/{bed}.bed",
#        bed_file = lambda wildcards: config["bedgz_annotation"][wildcards.sample] # can either be bed or bedgz
#    params:
#        chunk_size = 500000,
#        len_coulmn = lambda wildcards: config["bedgz_length_column"][wildcards.sample],
#    output:
#        flen_count_in_moitf = "flen_count_matrix/{sample}_intersect_{bed}_flen_in_motif.tsv",
#        flen_count_in_moitf_verbose = "flen_count_matrix/{sample}_intersect_{bed}_flen_in_motif.tsvverbos.tsv.gz"
#    shell:
#        "Rscript scripts/cfdna_center_in_tfbs.R {input.tfbs} {input.bed_file} {output.flen_count_in_moitf} {params.chunk_size} {params.len_coulmn} 1"    
        
rule count_flen_in_tfbs_new: 
    input:
        tfbs = "input_bed/{bed}.bed", 
        bed_file = lambda wildcards: config["bedgz_annotation"][wildcards.sample] # Make sure that cfDNA bed file is sorted
    params:
    output:
        flen_count_in_motif = "flen_count_matrix/{sample}_intersect_{bed}_flen_in_motif.tsv", 
        flen_count_in_motif_verbose = "flen_count_matrix/{sample}_intersect_{bed}_flen_in_motif.tsvverbos.tsv.gz", 
        center_map_on_bed = "flen_count_matrix/{sample}_intersect_{bed}_center_map.bed.gz", 

    shell:
        "sh scripts/map_cfdna_centers_to_tfbs.sh {input.tfbs} {input.bed_file}"
        " {output.flen_count_in_motif} {output.flen_count_in_motif_verbose}"
        " {output.center_map_on_bed}" 

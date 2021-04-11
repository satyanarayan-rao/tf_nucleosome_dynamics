rule generate_starting_head_tail_bed:
    input:
        head_tsv = "er_chia_pet_data/{sample}_head_location.tsv",
        tail_tsv = "er_chia_pet_data/{sample}_tail_location.tsv",
    params:
    output:
        starting_hg18_bed = "er_chia_pet_data/{sample}_head_tail_hg18.bed"
    shell:
        "sh scripts/chia_pet_prepare_bed.sh {input.head_tsv} {input.tail_tsv}"
        " {output.starting_hg18_bed}"
rule prepare_peak_location_bed:
    input:
        peak_location_tsv = "er_chia_pet_data/{sample}_peak_location.tsv"
    params:
    output:
        peak_location_bed = "er_chia_pet_data/{sample}_peak_location_hg18.bed"
    shell:
        "sh scripts/prepare_chia_pet_peak_bed.sh {input.peak_location_tsv}"
        " {output.peak_location_bed}"

rule liftover_hg18_to_hg38:
    input:
        starting_hg18_bed = "er_chia_pet_data/{sample}_head_tail_hg18.bed", 
        liftover_file = "/beevol/home/satyanarr/data/ucsc/liftover_chain_files/hg18ToHg38.over.chain"
    params:
    output:
        starting_hg38_bed = "er_chia_pet_data/{sample}_head_tail_hg38.bed"
    shell:
        "liftOver {input.starting_hg18_bed} {input.liftover_file} {output.starting_hg38_bed}"
        " \"{output.starting_hg38_bed}.tmp\"; rm \"{output.starting_hg38_bed}.tmp\" " 

rule combine_two_samples:
    input: 
        starting_hg38_bed1 = "er_chia_pet_data/{sample1}_head_tail_hg38.bed",
        starting_hg38_bed2 = "er_chia_pet_data/{sample2}_head_tail_hg38.bed"
    params:
    output:
        combined = "er_chia_pet_data/{sample1}_and_{sample2}_head_tail_combined_hg38.bed"
    shell:
        "cat {input.starting_hg38_bed1} {input.starting_hg38_bed2} > {output.combined}"

rule intersect_with_cnr:
    input:
        combined = "er_chia_pet_data/{sample1}_and_{sample2}_head_tail_combined_hg38.bed", 
        cnr_bed = lambda wildcards: config["cutnrun_data"][wildcards.cnr_bed]
    params:
    output: 
        intersect_with_cnr = "er_chia_pet_data/intersect_with_cnr/{sample1}_and_{sample2}_on_{cnr_bed}.bed"
    shell:
        "bedtools intersect -a {input.cnr_bed} -b {input.combined} -wa -wb > "
        " {output.intersect_with_cnr}"


rule extract_the_other_end_of_chia_pet:
    input:
        intersect_with_cnr = "er_chia_pet_data/intersect_with_cnr/{sample1}_and_{sample2}_on_{cnr_bed}.bed", 
        liftover_file = "/beevol/home/satyanarr/data/ucsc/liftover_chain_files/hg18ToHg38.over.chain"
    params:
    output:
        other_end = "er_chia_pet_data/other_end/{sample1}_and_{sample2}_on_{cnr_bed}_other_end.bed"
    shell:
        "sh scripts/extract_and_liftover_head_tail.sh {input.intersect_with_cnr}"
        " {input.liftover_file} {output.other_end}"

rule intersect_other_end_with_womotif:
    input:
        other_end = "er_chia_pet_data/other_end/{sample1}_and_{sample2}_on_{cnr_bed}_other_end.bed", 
        womotif = lambda wildcards: config["cutnrun_data"][wildcards.womotif]
    params:
    output:
        intersect = "er_chia_pet_data/other_end/{sample1}_and_{sample2}_on_{cnr_bed}_intersect_{womotif}.bed"
    shell:
        "bedtools intersect -a {input.womotif} -b {input.other_end} -wa -wb >"
        " {output.intersect}"
   

rule select_for_equal_number_of_random_class:
    input:
        observed = lambda wildcards: config["cutnrun_data"][wildcards.cls], 
        random_set = lambda wildcards: config["random_set"][wildcards.cls] 
    output:
        random_cls = "er_chia_pet_data/random_sites/sites_{cls}.bed"
    shell:
        "sh scripts/choose_randomly.sh {input.observed} {input.random_set} {output.random_cls}" 



#rule   
        
    
#rule get_corre


#### Archive ####### 

rule intersect_head_tail_to_peak:
    input: 
        starting_hg18_bed = "er_chia_pet_data/{sample}_head_tail_hg18.bed", 
        peak_location_bed = "er_chia_pet_data/{sample}_peak_location_hg18.bed"
    params:
    output: 
        peak_intersected_to_head_tail = "er_chia_pet_data/binding_sites_in_head_tail/{sample}_with_binding_sites_head_tail.bed"
    shell:
        "bedtools intersect -a {input.starting_hg18_bed} -b {input.peak_location_bed}"
        " -wa -wb > {output.peak_intersected_to_head_tail}"
    


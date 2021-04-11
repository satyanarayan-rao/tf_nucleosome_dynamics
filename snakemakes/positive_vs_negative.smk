rule fragment_heatmap_for_pos_and_neg:
    input:
        bed_file = "input_bed/{bed}.bed", 
        bw_file = lambda wildcards: config["chip_bigwigs"][wildcards.chip]
    params:
    output:
        
    shell:












rule fragment_heatmap_for_pos_and_neg:
    input:
    params:
    output:
    shell:


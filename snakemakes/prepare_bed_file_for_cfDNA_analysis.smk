rule prepare_bed_for_flank:
    input:
        start_bed = lambda wildcards: config[wildcards.bed_annotation]["start_bed"]
    params:
        annotation = lambda wildcards: config[wildcards.bed_annotation]["annotation"]
    output:
        bed_file = "input_bed/{bed_annotation}/{bed_annotation}_start.bed" 
    shell:
        "sh scripts/prepare_start_bed.sh {input.start_bed} {output.bed_file}"
        " {wildcards.bed_annotation} {params.annotation}"

rule flank_start_bed:
    input:
        bed_file = "input_bed/{bed_annotation}/{bed_annotation}_start.bed", 
        genome_size_file = lambda wildcards: config[wildcards.bed_annotation]["genome_size"], 
    params:
    output:
        slop_bed = "input_bed/{bed_annotation}/{bed_annotation}_slop_{slop}bp.bed", 
    shell:
        "sh scripts/slop_and_organize_name.sh {input.bed_file} {input.genome_size_file}"
        " {wildcards.slop} {output.slop_bed}"

rule flank_start_fasta:
    input:
        slop_bed = "input_bed/{bed_annotation}/{bed_annotation}_slop_{slop}bp.bed", 
        genome_fasta_file = lambda wildcards: config[wildcards.bed_annotation]["genome_fasta"], 
        
    params:
        bedtools_params = "-name"
    output:
        slop_fasta = "input_bed/{bed_annotation}/{bed_annotation}_slop_{slop}bp.fasta", 
    shell:
        "bedtools getfasta -fi {input.genome_fasta_file} -fo {output.slop_fasta}"
        " -bed {input.slop_bed} {params.bedtools_params}"


rule run_fimo:
    input:
        slop_fasta = "input_bed/{bed_annotation}/{bed_annotation}_slop_{slop}bp.fasta",
        meme_file = lambda wildcards: config[wildcards.bed_annotation]["motif_meme_file"] 
    params:
        fimo_params = "--max-strand --max-stored-scores 10000000" ,
        dir_name = lambda wildcards: "input_bed/{b}/fimo_on_{b}_slop_{s}bp".format(
                  b = wildcards.bed_annotation, s = wildcards.slop),
    output:
        fimo_out = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo.tsv"  
    shell:
        "fimo {params.fimo_params} --oc {params.dir_name}"
        " {input.meme_file} {input.slop_fasta}"


rule fimo2bed:
    input:
        fimo_out = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo.tsv", 
    params:
    output:
        fimo_bed = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo.bed", 
    shell:
        "sh scripts/process_fimo_and_convert_to_bed.sh {input.fimo_out}"
        " {output.fimo_bed}" 

# now slop by 50 bp

rule flank_fimo_bed:
    input:    
        fimo_bed = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo.bed", 
        genome_size_file = lambda wildcards: config[wildcards.bed_annotation]["genome_size"], 
    output:
        fimo_flanked = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_flank_{flank}bp.bed"
    shell:
        "bedtools slop -b {wildcards.flank} -g {input.genome_size_file}"
        " -i {input.fimo_bed} > {output.fimo_flanked}" 

rule remove_overlapping_regions_keep_best_by_score:
    input:
        fimo_flanked = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_flank_{flank}bp.bed"
    params:
        # from motif center keep the best motif in +-50 if overlaps
        non_overlapping_window = lambda wildcards: max(0, 100 - 2*int(wildcards.flank))
    output:
        non_overlapping_removed_best_kept = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/nonoverlapping_flank_{flank}bp.bed", 
        overlapping = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/overlapping_flank_{flank}bp.bed", 

    shell:
        "python scripts/prepare_non_overlapping_tfbs_bed.py {input.fimo_flanked}"
        " {params.non_overlapping_window} {output.non_overlapping_removed_best_kept}"
        " > {output.overlapping}" 
rule remove_chip_derived_blacklisted_regions:
    input:
        non_overlapping_removed_best_kept = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/nonoverlapping_flank_{flank}bp.bed",  
        blacklisted_region = lambda wildcards: config[wildcards.bed_annotation]["blacklisted_regions"]
    params:
    output:
        blacklisted_removed_best_kept = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/blacklist_removed_flank_{flank}bp.bed",  
        blacklisted = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/blacklisted_flank_{flank}bp.bed",  
        
    shell:
       "bedtools intersect -a {input.non_overlapping_removed_best_kept}"
       " -b {input.blacklisted_region} -wa -v > {output.blacklisted_removed_best_kept}"
       "; bedtools intersect -a {input.non_overlapping_removed_best_kept}"
       " -b {input.blacklisted_region} -wa -wb > {output.blacklisted}"

rule remove_chrY_for_wMotif:
    input:
        non_overlapping_removed_best_kept = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/blacklist_removed_flank_{flank}bp.bed",  
    params:
    output:
        chrY_removed = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/chrY_removed_flank_{flank}bp.bed"
    shell:
        "grep -v \"^chrY\" {input.non_overlapping_removed_best_kept} > {output.chrY_removed} "
rule make_a_copy_with_convenient_name_wMotif:
    input:
        chrY_removed = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/chrY_removed_flank_{flank}bp.bed"
    output:
        final_name = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/{flank}bp_smk_wMotif_{bed_annotation}.bed"
    shell:
        "cp {input.chrY_removed} {output.final_name}"

rule generate_count_table_wMotif:
    input:
        start_bed = lambda wildcards: config[wildcards.bed_annotation]["start_bed"], 
        bed_file = "input_bed/{bed_annotation}/{bed_annotation}_start.bed",
        slop_bed = "input_bed/{bed_annotation}/{bed_annotation}_slop_{slop}bp.bed", 
        fimo_bed = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo.bed", 
        non_overlapping_removed_best_kept = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/nonoverlapping_flank_{flank}bp.bed", 
        blacklisted_removed_best_kept = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/blacklist_removed_flank_{flank}bp.bed",  
        chrY_removed = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/chrY_removed_flank_{flank}bp.bed"
    params:
    output:
        count_table = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/count_table_{flank}bp.tsv"
    shell:
        "sh scripts/count_lines_on_list_of_files.sh"
        " {input.start_bed} {input.bed_file} {input.slop_bed} {input.fimo_bed}"
        " {input.non_overlapping_removed_best_kept} {input.blacklisted_removed_best_kept}"
        " {input.chrY_removed} {output.count_table}" 

rule motif_to_peak_center_distance_distribution:
    input:
        final_name = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/{flank}bp_smk_wMotif_{bed_annotation}.bed", 
    params: 
    output:
        distance_file = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/{flank}bp_motif_to_peak_distance.tsv"
    shell:
        "python scripts/motif_to_peak_center_dist.py {input.final_name}"
        " {output.distance_file} {wildcards.slop} "
    

rule generate_distance_histogram_wMotif:
    input:
        distance_file = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/{flank}bp_motif_to_peak_distance.tsv"
    params:
        target = lambda wildcards: config[wildcards.bed_annotation]["target"]
    output:
        motif_to_peak_dist_pdf = "input_bed/{bed_annotation}/plots/wMotif_to_peak_dist_slop_{slop}_flank_{flank}.pdf", 
        motif_to_peak_dist_png = "input_bed/{bed_annotation}/plots/wMotif_to_peak_dist_slop_{slop}_flank_{flank}.png", 
        motif_to_peak_dist_strand_pdf = "input_bed/{bed_annotation}/plots/strand_wMotif_to_peak_dist_slop_{slop}_flank_{flank}.pdf", 
        motif_to_peak_dist_strand_png = "input_bed/{bed_annotation}/plots/strand_wMotif_to_peak_dist_slop_{slop}_flank_{flank}.png", 
        
    shell:
        "Rscript scripts/plot_motif_to_peak_center_dist.R {input.distance_file}"
        " {output.motif_to_peak_dist_pdf} {output.motif_to_peak_dist_png}"
        " {params.target} \"wMotif\" {wildcards.slop}"
        " {output.motif_to_peak_dist_strand_pdf} {output.motif_to_peak_dist_strand_png}"
    

############### woMotif analysis ############# 

rule prepare_bed_for_wo_motif_analsyis:
    input:
        fimo_flanked = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_flank_{flank}bp.bed", 
        slop_bed = "input_bed/{bed_annotation}/{bed_annotation}_slop_{slop}bp.bed", 
    params:
    output:
        peaks_wo_motifs = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_flank_{flank}bp_removed_from_{bed_annotation}.bed" 
    shell:
        "python scripts/subtract_peaks_with_motif.py {input.fimo_flanked}"
        " {input.slop_bed} {output.peaks_wo_motifs}"  

rule get_fasta_for_woMotif:
    input:
        peaks_wo_motifs = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_flank_{flank}bp_removed_from_{bed_annotation}.bed", 
        genome_fasta_file = lambda wildcards: config[wildcards.bed_annotation]["genome_fasta"], 
        
    params:
        bedtools_params = "-name"
    output:
        slop_fasta = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_flank_{flank}bp_removed_from_{bed_annotation}.fasta", 
    shell:
        "bedtools getfasta -fi {input.genome_fasta_file} -fo {output.slop_fasta}"
        " -bed {input.peaks_wo_motifs} {params.bedtools_params}"


rule run_fimo_on_wo_motif_with_all_memes: 
    input:
        slop_fasta = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_flank_{flank}bp_removed_from_{bed_annotation}.fasta", 
        meme_file = lambda wildcards: config[wildcards.bed_annotation]["all_tf_memes"] 
    params:
        fimo_params = "--max-strand --max-stored-scores 10000000" ,
        dir_name = lambda wildcards: "input_bed/{b}/fimo_on_{b}_slop_{s}bp/fimo_on_woMotif_flank_{flank}".format(
                  b = wildcards.bed_annotation, s = wildcards.slop, flank = wildcards.flank),
    output:
        fimo_out = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_on_woMotif_flank_{flank}/fimo.tsv"  
    shell:
        "fimo {params.fimo_params} --oc {params.dir_name}"
        " {input.meme_file} {input.slop_fasta}"

rule choose_best_motif_for_each_peak:
    input:
        fimo_out = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_on_woMotif_flank_{flank}/fimo.tsv"  
    params: 
        annotation = lambda wildcards: config[wildcards.bed_annotation]["annotation"]
        
    output:
        peaks_with_top_motif = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_on_woMotif_flank_{flank}/fimo_with_top_motif.bed"
    shell:
        "sh scripts/choose_top_motif.sh  {input.fimo_out} {output.peaks_with_top_motif}"
        " {params.annotation}"


rule flank_woMotif_to_same_width_as_wMotif:
    input:
        chrY_removed = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/chrY_removed_flank_{flank}bp.bed", # to find width
        peaks_with_top_motif = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_on_woMotif_flank_{flank}/fimo_with_top_motif.bed", 
        genome_size_file = lambda wildcards: config[wildcards.bed_annotation]["genome_size"], 
 
    params:
    output:
        same_width_woMotif = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_on_woMotif_flank_{flank}/fimo_with_top_motif_same_width.bed"
    shell:
        "sh scripts/make_same_width.sh {input.chrY_removed} {input.peaks_with_top_motif}"
        " {output.same_width_woMotif} {input.genome_size_file}"

rule nonoverlapping_woMotif:
    input:
        same_width_woMotif = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_on_woMotif_flank_{flank}/fimo_with_top_motif_same_width.bed"
    params:
        # from motif center keep the best motif in +-50 if overlaps
        non_overlapping_window = lambda wildcards: max (0, 100 - 2*int(wildcards.flank)) 
    output:
        non_overlapping_removed_best_kept = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_on_woMotif_flank_{flank}/nonoverlapping_woMotif_flank_{flank}bp.bed", 
        overlapping_removed_best_kept = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_on_woMotif_flank_{flank}/overlapping_woMotif_flank_{flank}bp.bed", 
    shell:
        "python scripts/prepare_non_overlapping_tfbs_bed.py {input.same_width_woMotif}"
        " {params.non_overlapping_window} {output.non_overlapping_removed_best_kept}"
        " > {output.overlapping_removed_best_kept}"

rule remove_chip_derived_blacklisted_regions_woMotif:
    input:
        non_overlapping_removed_best_kept = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_on_woMotif_flank_{flank}/nonoverlapping_woMotif_flank_{flank}bp.bed", 
        blacklisted_region = lambda wildcards: config[wildcards.bed_annotation]["blacklisted_regions"]
    params:
    output:
        blacklisted_removed_best_kept = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_on_woMotif_flank_{flank}/blacklist_removed_flank_{flank}bp.bed",  
        blacklisted = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_on_woMotif_flank_{flank}/blacklisted_flank_{flank}bp.bed",  
        
    shell:
       "bedtools intersect -a {input.non_overlapping_removed_best_kept}"
       " -b {input.blacklisted_region} -wa -v > {output.blacklisted_removed_best_kept}"
       "; bedtools intersect -a {input.non_overlapping_removed_best_kept}"
       " -b {input.blacklisted_region} -wa -wb > {output.blacklisted}"


rule remove_chrY_woMotif:
    input:
        non_overlapping_removed_best_kept = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_on_woMotif_flank_{flank}/blacklist_removed_flank_{flank}bp.bed", 
        
    params:
    output: 
        chrY_removed = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_on_woMotif_flank_{flank}/chrY_removed_flank_{flank}bp.bed", 
        
    shell:
        "grep -v \"^chrY\" {input.non_overlapping_removed_best_kept} > {output.chrY_removed}"


rule make_a_copy_with_convenient_name_woMotif:
    input:
        chrY_removed = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_on_woMotif_flank_{flank}/chrY_removed_flank_{flank}bp.bed", 
    output:
        final_name = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_on_woMotif_flank_{flank}/{flank}bp_smk_woMotif_{bed_annotation}.bed"
    shell:
        "cp {input.chrY_removed} {output.final_name}"

rule generate_count_table_woMotif:
    input:
        slop_bed = "input_bed/{bed_annotation}/{bed_annotation}_slop_{slop}bp.bed", 
        fimo_bed = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo.bed", 
        peaks_wo_motifs = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_flank_{flank}bp_removed_from_{bed_annotation}.bed",
        non_overlapping_removed_best_kept = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_on_woMotif_flank_{flank}/nonoverlapping_woMotif_flank_{flank}bp.bed", 
        blacklisted_removed_best_kept = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_on_woMotif_flank_{flank}/blacklist_removed_flank_{flank}bp.bed",  
        
        chrY_removed = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_on_woMotif_flank_{flank}/chrY_removed_flank_{flank}bp.bed", 
    params:
    output:
        count_table = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_on_woMotif_flank_{flank}/count_table_{flank}bp.tsv"
    shell:
        "sh scripts/count_lines_on_list_of_files.sh"
        " {input.slop_bed} {input.fimo_bed} {input.peaks_wo_motifs}"
        " {input.non_overlapping_removed_best_kept} {input.blacklisted_removed_best_kept}"
        " {input.chrY_removed} {output.count_table}"
        "; sh scripts/get_unique_peaks_from_fimo.sh {input.fimo_bed} {output.count_table}"

rule motif_to_peak_center_distance_distribution_woMotif:
    input:
        final_name = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_on_woMotif_flank_{flank}/{flank}bp_smk_woMotif_{bed_annotation}.bed"
    params: 
    output:
        distance_file = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_on_woMotif_flank_{flank}/woMotif_motif_to_peak_distance.tsv"
    shell:
        "python scripts/motif_to_peak_center_dist.py {input.final_name}"
        " {output.distance_file} {wildcards.slop}"

rule generate_distance_histogram_woMotif:
    input:
        distance_file = "input_bed/{bed_annotation}/fimo_on_{bed_annotation}_slop_{slop}bp/fimo_on_woMotif_flank_{flank}/woMotif_motif_to_peak_distance.tsv"
    params:
        target = lambda wildcards: config[wildcards.bed_annotation]["target"]
    output:
        motif_to_peak_dist_pdf = "input_bed/{bed_annotation}/plots/woMotif_to_peak_dist_slop_{slop}_flank_{flank}.pdf", 
        motif_to_peak_dist_png = "input_bed/{bed_annotation}/plots/woMotif_to_peak_dist_slop_{slop}_flank_{flank}.png",
        motif_to_peak_dist_strand_pdf = "input_bed/{bed_annotation}/plots/strand_woMotif_to_peak_dist_slop_{slop}_flank_{flank}.pdf",
        motif_to_peak_dist_strand_png = "input_bed/{bed_annotation}/plots/strand_woMotif_to_peak_dist_slop_{slop}_flank_{flank}.png",
    shell:
        "Rscript scripts/plot_motif_to_peak_center_dist.R {input.distance_file}"
        " {output.motif_to_peak_dist_pdf} {output.motif_to_peak_dist_png}"
        " {params.target} \"woMotif\" {wildcards.slop}"
        " {output.motif_to_peak_dist_strand_pdf} {output.motif_to_peak_dist_strand_png}" 

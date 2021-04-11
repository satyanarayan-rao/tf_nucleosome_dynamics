def  get_genome_size_file(wildcards):
    tfbs_bed_without_flank = "".join(wildcards.bed.split("bp_")[1:])
    tfbs_bed_without_flank = tfbs_bed_without_flank.replace(".bed", "") 
    return config["tfbs_genome_annotation"][tfbs_bed_without_flank] 

rule generate_tfbs_flank_for_vplot: 
    input:
        tfbs = "input_bed/{bed}.bed",
        genome_size_file = lambda wildcards: get_genome_size_file(wildcards)
    params: 
        to_slop = lambda wildcards: int(wildcards.flank) - int (wildcards.bed.split("bp_")[0])
    output:
        vplot_bed = "vplot_input_bed/{bed}_total_flank_{flank}.bed"
    shell:
        "bedtools slop -b {params.to_slop} -g {input.genome_size_file}"
        " -i {input.tfbs} > {output.vplot_bed}" 
#rule get_data_for_vplot_matrix_generation:
#    input:
#        vplot_bed = "vplot_input_bed/{bed}_total_flank_{flank}.bed", 
#        bed_file = lambda wildcards: config["bedgz_annotation"][wildcards.sample] # can either be bed or bedgz
#    params:
#        chunk_size = 500000,
#        len_coulmn = lambda wildcards: config["bedgz_length_column"][wildcards.sample],
#    output:
#        flen_count_in_moitf = "vplot_flen_count_matrix/{sample}_intersect_{bed}_total_flank_{flank}_flen_in_motif.tsv", 
#        flen_count_in_moitf_verbose = "vplot_flen_count_matrix/{sample}_intersect_{bed}_total_flank_{flank}_flen_in_motif.tsvverbos.tsv.gz"
#    shell:
#        "Rscript scripts/cfdna_center_in_tfbs.R {input.vplot_bed} {input.bed_file} {output.flen_count_in_moitf} {params.chunk_size} {params.len_coulmn} 1"  

## Implementing new rule that uses bedtools to intersect

rule  get_data_for_vplot_matrix_generation:
    input:
        vplot_bed = "vplot_input_bed/{bed}_total_flank_{flank}.bed", 
        bed_file = lambda wildcards: config["bedgz_annotation"][wildcards.sample] # can either be bed or bedgz
    params:
    output:
        flen_count_in_motif = "vplot_flen_count_matrix/{sample}_intersect_{bed}_total_flank_{flank}_flen_in_motif.tsv", 
        flen_count_in_motif_verbose = "vplot_flen_count_matrix/{sample}_intersect_{bed}_total_flank_{flank}_flen_in_motif.tsvverbos.tsv.gz", 
        center_map_on_bed = "vplot_flen_count_matrix/{sample}_intersect_{bed}_total_flank_{flank}_center_map.bed.gz"
    shell:
        "sh scripts/map_cfdna_centers_to_tfbs.sh {input.vplot_bed} {input.bed_file}"
        " {output.flen_count_in_motif} {output.flen_count_in_motif_verbose}"
        " {output.center_map_on_bed}"

rule gen_vplot_matrix_for_all_sites:
    input:
        flen_count_in_motitf_verbose_from_vplot = "vplot_flen_count_matrix/{sample}_intersect_{bed}_total_flank_{flank}_flen_in_motif.tsvverbos.tsv.gz", 
        flen_count_in_motif_tfbs_flank = "flen_count_matrix/{sample}_intersect_{bed}_flen_in_motif.tsv",  # this will give me a freedom to select based on minimum number of fragments centers I have 
        start_bed_file = "input_bed/{bed}.bed"
    params:
        len_column = lambda wildcards: config["bedgz_length_column"][wildcards.sample], 
        actual_flank = lambda wildcards: wildcards.bed.split("bp_")[0]
    output:
        vplot_matrix_based_on_count_in_tfbs = "vplot_matrix_for_sites_with_th_count/{sample}_intersect_{bed}_total_flank_{flank}_tfbs_th_count_{cnt}_in_range_{min_frag_len}_to_{max_frag_len}.tsv", 
        tfbs_passing_th = "sites_passing_dna_count_th/{sample}_intersect_{bed}_total_flank_{flank}_tfbs_th_count_{cnt}_in_range_{min_frag_len}_to_{max_frag_len}.tsv"
        
    shell:
        "sh scripts/gen_vplot_data_tfbs_th.sh {input.flen_count_in_motitf_verbose_from_vplot}" 
        " {output.vplot_matrix_based_on_count_in_tfbs} {input.start_bed_file}"
        " {params.len_column} {input.flen_count_in_motif_tfbs_flank}"
        " {wildcards.cnt} {wildcards.flank} {params.actual_flank}" 
        " {wildcards.min_frag_len} {wildcards.max_frag_len}"
        " {output.tfbs_passing_th}"  

rule generate_vplot_for_all_sites:
    input:
        vplot_matrix_based_on_count_in_tfbs = "vplot_matrix_for_sites_with_th_count/{sample}_intersect_{bed}_total_flank_{flank}_tfbs_th_count_{cnt}_in_range_{min_frag_len}_to_{max_frag_len}.tsv",  
        tfbs_passing_th = "sites_passing_dna_count_th/{sample}_intersect_{bed}_total_flank_{flank}_tfbs_th_count_{cnt}_in_range_{min_frag_len}_to_{max_frag_len}.tsv", 
        start_tfbs = "input_bed/{bed}.bed"
         
    params:
        bed_file_annotation = lambda wildcards: config["bed_annotation"][wildcards.bed.split("bp_")[1].replace(".bed","")], 
        cbra_lim = lambda wildcards: config["vplot_cbra_limits"][wildcards.sample]
    output:
        vplot_pdf = "plots/vplot_for_all_sites/{sample}_intersect_{bed}_total_flank_{flank}_tfbs_th_count_{cnt}_in_range_{min_frag_len}_to_{max_frag_len}.pdf", 
        vplot_gplt = "plots/vplot_for_all_sites/{sample}_intersect_{bed}_total_flank_{flank}_tfbs_th_count_{cnt}_in_range_{min_frag_len}_to_{max_frag_len}.gplt",
        vplot_zoom_pdf = "plots/vplot_for_all_sites/{sample}_intersect_{bed}_total_flank_{flank}_tfbs_th_count_{cnt}_in_range_{min_frag_len}_to_{max_frag_len}_vplot_zoomed.pdf", 
        vplot_zoom_gplt = "plots/vplot_for_all_sites/{sample}_intersect_{bed}_total_flank_{flank}_tfbs_th_count_{cnt}_in_range_{min_frag_len}_to_{max_frag_len}_vplot_zoomed.gplt", 
        layered_zoom_pdf = "plots/vplot_for_all_sites/{sample}_intersect_{bed}_total_flank_{flank}_tfbs_th_count_{cnt}_in_range_{min_frag_len}_to_{max_frag_len}_vplot_zoomed_split.pdf",  
        layered_zoom_gplt = "plots/vplot_for_all_sites/{sample}_intersect_{bed}_total_flank_{flank}_tfbs_th_count_{cnt}_in_range_{min_frag_len}_to_{max_frag_len}_vplot_zoomed_split.gplt",
        zscore_tsv = "plots/vplot_for_all_sites/{sample}_intersect_{bed}_total_flank_{flank}_tfbs_th_count_{cnt}_in_range_{min_frag_len}_to_{max_frag_len}_vplot_zoomed_zscore.tsv",
        zscored_zoom_gplt = "plots/vplot_for_all_sites/{sample}_intersect_{bed}_total_flank_{flank}_tfbs_th_count_{cnt}_in_range_{min_frag_len}_to_{max_frag_len}_vplot_zoomed_zscore.gplt", 
        zscored_zoom_whole_pdf = "plots/vplot_for_all_sites/{sample}_intersect_{bed}_total_flank_{flank}_tfbs_th_count_{cnt}_in_range_{min_frag_len}_to_{max_frag_len}_vplot_zoomed_zscore.pdf",
        zscored_zoom_splt_pdf = "plots/vplot_for_all_sites/{sample}_intersect_{bed}_total_flank_{flank}_tfbs_th_count_{cnt}_in_range_{min_frag_len}_to_{max_frag_len}_vplot_zoomed_zscore_split.pdf", 
    shell:
        "sh scripts/gen_vplot_all_sites.sh {input.vplot_matrix_based_on_count_in_tfbs}"
        " {output.vplot_gplt} {output.vplot_pdf} {wildcards.flank}"
        " {wildcards.sample} {params.bed_file_annotation} \"all_sites\"" 
        " {output.vplot_zoom_gplt} {output.vplot_zoom_pdf} {params.cbra_lim}"
        " {input.tfbs_passing_th} {output.layered_zoom_gplt} {output.layered_zoom_pdf}" 
        " {output.zscore_tsv} {output.zscored_zoom_gplt} {output.zscored_zoom_whole_pdf}" 
        " {output.zscored_zoom_splt_pdf} {input.start_tfbs}"
rule assing_cluster_id_to_vplot_data:
    input:
        flen_count_in_moitf_verbose = "vplot_flen_count_matrix/{sample}_intersect_{bed}_total_flank_{flank}_flen_in_motif.tsvverbos.tsv.gz", 
        reassigned_bed_tsv = "fragment_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_reassigned_bed.tsv",
    params: 
        actual_flank = lambda wildcards: wildcards.bed.split("bp_")[0]
    output:
        flen_count_in_motif_with_cl_id = "cluster_id_tagged_to_vplot_flen_count_matrix/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_flen.tsv.gz"
    shell: 
        "python scripts/tag_cl_id_for_vplot.py {input.flen_count_in_moitf_verbose}"
        " {input.reassigned_bed_tsv} {output.flen_count_in_motif_with_cl_id} {wildcards.flank} {params.actual_flank}" 

rule generate_vplot_matrix_for_a_cluster:
    input:
        flen_count_in_motif_with_cl_id = "cluster_id_tagged_to_vplot_flen_count_matrix/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_flen.tsv.gz", 
        start_bed_file = "input_bed/{bed}.bed", 
    params: 
        len_coulmn = lambda wildcards: config["bedgz_length_column"][wildcards.sample],
    output:
        vplot_matrix_for_cfdna_cl = "cluster_specific_vplot_matrix/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot.tsv", 
    shell:
        "sh scripts/get_cluster_wise_vplot_data.sh {input.flen_count_in_motif_with_cl_id}"
        " {wildcards.cl_id} {output.vplot_matrix_for_cfdna_cl}"
        " {input.start_bed_file} {params.len_coulmn}"

rule prepend_zeors_to_vplot_matrix:
    input:
        vplot_matrix_for_cfdna_cl = "cluster_specific_vplot_matrix/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot.tsv", 
       
    params:
    output:
        zero_prepended_matrix  = "cluster_specific_vplot_matrix/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot_zero_pad.tsv", 
    shell:
        "python scripts/pad_zeros_to_vplot_matrix.py {input.vplot_matrix_for_cfdna_cl} {output.zero_prepended_matrix}"


rule generate_vplot:
    input:
        #vplot_matrix_for_cfdna_cl = "cluster_specific_vplot_matrix/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot.tsv",
        zero_prepended_matrix  = "cluster_specific_vplot_matrix/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot_zero_pad.tsv",  
        reassigned_bed_tsv = "fragment_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_reassigned_bed.tsv",
        
    params:
        bed_file_annotation = lambda wildcards: config["bed_annotation"][wildcards.bed.split("bp_")[1].replace(".bed","")], 
        cbra_lim = lambda wildcards: config["vplot_cbra_limits"][wildcards.sample]
    output: 
        vplot_cl_pdf = "plots/vplot_for_cfdna_cl/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot.pdf",
        vplot_cl_gplt = "plots/vplot_for_cfdna_cl/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot.gplt", 
        vplot_cl_zoom_pdf = "plots/vplot_for_cfdna_cl/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot_zoomed.pdf", 
        vplot_cl_zoom_gplt = "plots/vplot_for_cfdna_cl/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot_zoomed.gplt", 
         layered_cl_zoom_pdf = "plots/vplot_for_cfdna_cl/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot_zoomed_split.pdf",  
         layered_cl_zoom_gplt = "plots/vplot_for_cfdna_cl/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot_zoomed_split.gplt",
         zscore_tsv = "plots/vplot_for_cfdna_cl/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot_zoomed_zscore.tsv",
         zscored_zoom_gplt = "plots/vplot_for_cfdna_cl/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot_zoomed_zscore.gplt", 
         zscored_zoom_whole_pdf = "plots/vplot_for_cfdna_cl/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot_zoomed_zscore.pdf",
         zscored_zoom_splt_pdf = "plots/vplot_for_cfdna_cl/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot_zoomed_zscore_split.pdf", 
 
        
    shell:
        #"sh scripts/gen_vplot.sh {input.vplot_matrix_for_cfdna_cl}"
        "sh scripts/gen_vplot.sh {input.zero_prepended_matrix}"
        " {output.vplot_cl_gplt} {output.vplot_cl_pdf} {wildcards.flank}"
        " {wildcards.sample} {params.bed_file_annotation} {wildcards.cl_id}"
        " {output.vplot_cl_zoom_gplt} {output.vplot_cl_zoom_pdf} {params.cbra_lim}"
        " {input.reassigned_bed_tsv} {output.layered_cl_zoom_gplt} {output.layered_cl_zoom_pdf}" 
        " {output.zscore_tsv} {output.zscored_zoom_gplt} {output.zscored_zoom_whole_pdf}" 
        " {output.zscored_zoom_splt_pdf}"

    

#### generate C&R vplots for length clusters #### 


## generate cnr v matrix data

rule get_cnr_data_for_vplot_matrix_generation:
    input:
        vplot_bed = "vplot_input_bed/{bed}_total_flank_{flank}.bed", 
        bed_file = lambda wildcards: config["cnr_bedgz_annotation"][wildcards.cnr] # can either be bed or bedgz
    params:
        chunk_size = 500000,
        len_coulmn = lambda wildcards: config["cnr_bedgz_length_column"][wildcards.cnr],
    output:
        flen_count_in_moitf = "cnr_vplot_flen_count_matrix/{cnr}_intersect_{bed}_total_flank_{flank}_flen_in_motif.tsv", 
        flen_count_in_moitf_verbose = "cnr_vplot_flen_count_matrix/{cnr}_intersect_{bed}_total_flank_{flank}_flen_in_motif.tsvverbos.tsv.gz"
    shell:
        "Rscript scripts/cfdna_center_in_tfbs.R {input.vplot_bed} {input.bed_file} {output.flen_count_in_moitf} {params.chunk_size} {params.len_coulmn} 1"  

rule assign_cluster_id_to_cnr_vplot_data:
    input:
        flen_count_in_moitf_verbose = "cnr_vplot_flen_count_matrix/{cnr}_intersect_{bed}_total_flank_{flank}_flen_in_motif.tsvverbos.tsv.gz", 
        reassigned_bed_tsv = "fragment_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_reassigned_bed.tsv",
    params:
        actual_flank = lambda wildcards: wildcards.bed.split("bp_")[0]
    output:
        flen_count_in_motif_with_cl_id = "cluster_id_tagged_to_cnr_vplot_flen_count_matrix/{cnr}_on_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_flen.tsv.gz"
    shell:
        "python scripts/tag_cl_id_for_vplot.py {input.flen_count_in_moitf_verbose}"
        " {input.reassigned_bed_tsv} {output.flen_count_in_motif_with_cl_id} {wildcards.flank} {params.actual_flank}" 

rule generate_cnr_vplot_matrix_for_a_cluster:
    input:
        flen_count_in_motif_with_cl_id = "cluster_id_tagged_to_cnr_vplot_flen_count_matrix/{cnr}_on_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_flen.tsv.gz", 
        start_bed_file = "input_bed/{bed}.bed", 
    params: 
        #len_coulmn = l # because for cnr data it is always 4 #lambda wildcards: config["bedgz_length_column"][wildcards.sample],
        len_coulmn = lambda wildcards: config["cnr_bedgz_length_column"][wildcards.cnr],
    output:
        vplot_matrix_for_cfdna_cl = "cluster_specific_cnr_vplot_matrix/{cnr}_on_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot.tsv", 
    shell:
        "sh scripts/get_cluster_wise_vplot_data.sh {input.flen_count_in_motif_with_cl_id}"
         " {wildcards.cl_id} {output.vplot_matrix_for_cfdna_cl}"
         " {input.start_bed_file} {params.len_coulmn}"

rule prepend_zeors_to_cnr_vplot_matrix:
    input:
        vplot_matrix_for_cfdna_cl = "cluster_specific_cnr_vplot_matrix/{cnr}_on_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot.tsv", 
       
    params:
    output:
        zero_prepended_matrix  = "cluster_specific_cnr_vplot_matrix/{cnr}_on_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot_zero_pad.tsv", 
    shell:
        "python scripts/pad_zeros_to_vplot_matrix.py {input.vplot_matrix_for_cfdna_cl} {output.zero_prepended_matrix}"

rule generate_cnr_vplot:
    input:
        #vplot_matrix_for_cfdna_cl = "cluster_specific_cnr_vplot_matrix/{cnr}_on_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot.tsv",
        zero_prepended_matrix  = "cluster_specific_cnr_vplot_matrix/{cnr}_on_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot_zero_pad.tsv", 
        reassigned_bed_tsv = "fragment_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_reassigned_bed.tsv",
    params:
        bed_file_annotation = lambda wildcards: config["bed_annotation"][wildcards.bed.split("bp_")[1].replace(".bed","")], 
        cbra_lim = lambda wildcards: config["cnr_vplot_cbra_limits"][wildcards.sample], 
         #layout_lims = lambda wildcards: "@".join(map(str, config["cnr_plot_arrangment"][wildcards.layout]))
    output: 
        vplot_cl_pdf = "plots/cnr_vplot_for_cfdna_cl/{cnr}_on_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot.pdf",
        vplot_cl_gplt = "plots/cnr_vplot_for_cfdna_cl/{cnr}_on_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot.gplt", 
        vplot_cl_zoom_pdf = "plots/cnr_vplot_for_cfdna_cl/{cnr}_on_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot_zoomed.pdf", 
        vplot_cl_zoom_gplt = "plots/cnr_vplot_for_cfdna_cl/{cnr}_on_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot_zoomed.gplt",  
         layered_cl_zoom_pdf = "plots/cnr_vplot_for_cfdna_cl/{cnr}_on_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot_zoomed_split.pdf",  
         layered_cl_zoom_gplt = "plots/cnr_vplot_for_cfdna_cl/{cnr}_on_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot_zoomed_split.gplt",
         zscore_tsv = "plots/cnr_vplot_for_cfdna_cl/{cnr}_on_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot_zoomed_zscore.tsv",
         zscored_zoom_gplt = "plots/cnr_vplot_for_cfdna_cl/{cnr}_on_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot_zoomed_zscore.gplt", 
         zscored_zoom_whole_pdf = "plots/cnr_vplot_for_cfdna_cl/{cnr}_on_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot_zoomed_zscore.pdf",
         zscored_zoom_splt_pdf = "plots/cnr_vplot_for_cfdna_cl/{cnr}_on_{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_flank_{flank}_cl_{cl_id}_vplot_zoomed_zscore_split.pdf", 
         
    shell:
        #"sh scripts/gen_cnr_vplot.sh {input.vplot_matrix_for_cfdna_cl}"
        "sh scripts/gen_cnr_vplot.sh {input.zero_prepended_matrix}"
        " {output.vplot_cl_gplt} {output.vplot_cl_pdf} {wildcards.flank}"
        " {wildcards.sample} {params.bed_file_annotation} {wildcards.cl_id}"
        " {output.vplot_cl_zoom_gplt} {output.vplot_cl_zoom_pdf} {params.cbra_lim}"
        " {input.reassigned_bed_tsv} {output.layered_cl_zoom_gplt} {output.layered_cl_zoom_pdf}" 
        " {output.zscore_tsv} {output.zscored_zoom_gplt} {output.zscored_zoom_whole_pdf}" 
        " {output.zscored_zoom_splt_pdf}"


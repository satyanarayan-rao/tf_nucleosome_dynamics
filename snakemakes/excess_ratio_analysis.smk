rule reformat_bed_to_tsv: 
    input:
        tsv_file = "fragment_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_reassigned_bed.tsv",
    output:
        bed_file = "fragment_bw_map_manual/{sample}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_reassigned.bed",
    shell:
        "sh scripts/tsv2bed.sh {input.tsv_file} {output.bed_file}"

rule matrix_for_sites_in_length_clusters:
    input:
        common_bed_a = "fragment_bw_map_manual/{sample_a}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_reassigned.bed",
        common_bed_b = "fragment_bw_map_manual/{sample_b}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_reassigned.bed",
    params:
    output:
        table_cl_ij = "tfbs_cnt_matrix_for_excess_ratio/{sample_a}_vs_{sample_b}_intersect_{bed}_dm_flen_min_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_nclust_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cl_{ref_for_ks_test}_table.tsv"
    shell:
        "sh scripts/gen_tfbs_count_matrix_a_vs_b.sh {input.common_bed_a}"
        " {input.common_bed_b} {output.table_cl_ij}" 
rule generate_excess_ratio_matrix: 
    input:
        table_cl_ij = "tfbs_cnt_matrix_for_excess_ratio/{sample_a}_vs_{sample_b}_intersect_{bed}_dm_flen_min_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_nclust_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cl_{ref_for_ks_test}_table.tsv"
    params:
    output:
        excess_ratio_matrix = "excess_ratio_matrix/{sample_a}_vs_{sample_b}_intersect_{bed}_dm_flen_min_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_nclust_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cl_{ref_for_ks_test}_excess_ratio.tsv", 
        adjusted_p_val_matrix = "excess_ratio_matrix/{sample_a}_vs_{sample_b}_intersect_{bed}_dm_flen_min_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_nclust_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cl_{ref_for_ks_test}_adj_p_val.tsv",
        adjusted_s_str_matrix = "excess_ratio_matrix/{sample_a}_vs_{sample_b}_intersect_{bed}_dm_flen_min_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_nclust_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cl_{ref_for_ks_test}_adj_s_str.tsv",  
    shell:
        "Rscript scripts/calculate_excess_ratio.R {input.table_cl_ij}"
        " {output.excess_ratio_matrix}"
        "; Rscript scripts/do_phyper_of_overlap.R"
        " {input.table_cl_ij} {output.adjusted_p_val_matrix}"
        " {output.adjusted_s_str_matrix} {output.excess_ratio_matrix}" 
    
rule plot_excess_ratio_matrix:
    input: 
        excess_ratio_matrix = "excess_ratio_matrix/{sample_a}_vs_{sample_b}_intersect_{bed}_dm_flen_min_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_nclust_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cl_{ref_for_ks_test}_excess_ratio.tsv", 
        adjusted_s_str_matrix = "excess_ratio_matrix/{sample_a}_vs_{sample_b}_intersect_{bed}_dm_flen_min_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_nclust_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cl_{ref_for_ks_test}_adj_s_str.tsv",   
        base_excess_ratio_gplt = "utils/gnuplot_base_files/base_excess_ratio.gplt" 
        
    params:
    output:
        excess_ratio_pdf = "plots/excess_ratio_heatmaps/{sample_a}_vs_{sample_b}_intersect_{bed}_dm_flen_min_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_nclust_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cl_{ref_for_ks_test}_excess_ratio.pdf",
        excess_ratio_gplt = "plots/excess_ratio_heatmaps/{sample_a}_vs_{sample_b}_intersect_{bed}_dm_flen_min_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_nclust_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cl_{ref_for_ks_test}_excess_ratio.gplt",
    shell:
        "python scripts/generate_gnuplt_file_for_excess_matrix.py"
        " {input.excess_ratio_matrix} {input.base_excess_ratio_gplt}" 
        " {output.excess_ratio_gplt} {wildcards.nclust} {output.excess_ratio_pdf}"
        " {wildcards.sample_a} {wildcards.sample_b} {input.adjusted_s_str_matrix}"
        " ; gnuplot {output.excess_ratio_gplt}"  

# Subtracting significant sites from one system to the other ###  
# For example to define MCF-7 CDX specific ER E2 sites - we remove
# sites that belong to length clusters with 
# significantly higher C&R scores (compared to the nuc cluster - cluster 6) in IH02
# NOTE: here only common sites should be used. 
def get_files_for_subtraction(wildcards, target):
    flist = []
    # get reassigned bed for the target
    sample_list = config["subtract_significant_sites"][wildcards.bed][target]["sample"]
    for s in sample_list:
        fname = "fragment_bw_map_manual/{s}_intersect_{bed}_dmatrix_flen_min_{min_flen}_max_{max_flen}_least_read_{min_reads}_bw_{bw}_nclust_{nclust}_max_iter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cluster_{ref_for_ks_test}_reassigned.bed".format(
            s=s, bed = wildcards.bed, min_flen=wildcards.min_flen, 
            max_flen=wildcards.max_flen, min_reads=wildcards.min_reads,
            bw=wildcards.bw, nclust=wildcards.nclust, max_iter=wildcards.max_iter,
            slop=wildcards.slop, fsize=wildcards.fsize, chip=wildcards.chip,
            ref_for_ks_test=wildcards.ref_for_ks_test)
        flist.append(fname)
    return flist

def get_cl_id_for_subtraction(wildcards, target):
    cl_str_list = [] 
    # get reassigned bed for the target
    sig_cl_list = config["subtract_significant_sites"][wildcards.bed][target]["significant_cl"]
    for cl_list in sig_cl_list: # check `subtract_significant_sites` in configs/config.yaml
        cl_str = "@".join(map(str, cl_list))  
        cl_str_list.append(cl_str) 
    final_str = " ".join(cl_str_list)
    
    return final_str


rule subtract_significant_sites: 
    input:
        common_bed_a = lambda wildcards: get_files_for_subtraction(wildcards, wildcards.setting_a), 
        common_bed_b = lambda wildcards: get_files_for_subtraction(wildcards, wildcards.setting_b), 
    params: 
        cluster_vec_str_a = lambda wildcards: get_cl_id_for_subtraction(wildcards, wildcards.setting_a), 
        cluster_vec_str_b = lambda wildcards: get_cl_id_for_subtraction(wildcards, wildcards.setting_b), 
    output:
        subtracted_bed = "significant_sites_subtracted/{setting_b}_subtracted_from_{setting_a}_for_{bed}_flen_min_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cl_{ref_for_ks_test}_reassigned.bed", 
        common_bed = "significant_sites_subtracted/{setting_b}_common_{setting_a}_for_{bed}_flen_min_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cl_{ref_for_ks_test}_reassigned.bed"
    shell:
        "python scripts/subtract_significant_sites.py"
        " \"{input.common_bed_a}\" \"{params.cluster_vec_str_a}\""
        " \"{input.common_bed_b}\" \"{params.cluster_vec_str_b}\""
        " {output.subtracted_bed} {output.common_bed}" 
    

#rule map_to_the_specific_sites: 
#    input:
#        subtracted_bed = "significant_sites_subtracted/{setting_b}_subtracted_from_{setting_a}_for_{bed}_flen_min_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cl_{ref_for_ks_test}_reassigned.bed",   
#        bw_file = lambda wildcards: 
#    params:
#    output:
#        bw_map = "bigwig_mapped_to_sig_sites_subtracted/{bigwig}_to_{setting_b}_subtracted_from_{setting_a}_for_{bed}_flen_min_{min_flen}_max_{max_flen}_lread_{min_reads}_bw_{bw}_ncl_{nclust}_miter_{max_iter}_slop_{slop}_fsize_{fsize}_chip_{chip}_ref_cl_{ref_for_ks_test}_reassigned.csv.gz"
#        
#    shell:
    


#rule remove_significan_and_specific_sites:
#    input: 
#        
#    params:
#    output:
#    shell:

# Mapping Transcription Factor-Nucleosome Dynamics from Plasma cfDNA

Scripts and snakemake pipeline for generating manuscript figures. An example to generate CTCF ChIP score boxplot is shown below: 

```
$ snakemake --snakefile discover_tf_binding_patterns.smk plots/chip_box_plots/ih02_intersect_50bp_ctcf_enc_carroll_etoh_non_overlapping_dmatrix_flen_min_35_max_250_least_read_5_bw_3_nclust_6_max_iter_250_slop_300_chip_GM_ENCFF578TBN_ref_cluster_6_bxplt.pdf --configfile configs/config.yaml --cluster-config configs/cluster.json --cluster "bsub -q {cluster.queue} -R {cluster.resources} -J {cluster.name} -o {cluster.output} -e {cluster.error}"  -j50 
```

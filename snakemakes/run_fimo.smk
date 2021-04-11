rule run_fimo_on_hg38:
    input:
        fa = "/beevol/home/satyanarr/data/ucsc/human_genome/hg38/hg38.fa",
        meme = "/beevol/home/satyanarr/cfDNA-TF-footprints/TFs/motifs/ER/versions.meme"
    params:
        fimo_params = "--max-strand --max-stored-scores 10000000" ,
        dir_name = "er_motifs_in_hg38"
    output:
        fimo_out = "er_motifs_in_hg38/fimo.tsv"
    shell:
        "fimo {params.fimo_params} --oc {params.dir_name} {input.meme} {input.fa}" 


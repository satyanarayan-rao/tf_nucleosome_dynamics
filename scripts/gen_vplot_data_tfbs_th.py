import sys
from collections import defaultdict
import re
import gzip
from collections import defaultdict
# zcat flen_count_matrix/ih02_intersect_50bp_gm12878_ENCFF018NNF_non_overlapping_flen_in_motif.tsvverbos.tsv.gz  | head -1
# chr1	805021	805128	chr1:804614-805250|gm12878_ENCFF018NNF`27010_CTGGGAGG	11.8171	-	chr1	804891	805225	IH02_00.pairs@123159	334	+

verbose_fp = gzip.open(sys.argv[1])
fragment_len_column = int(sys.argv[2]) 
bed_cols = int(sys.argv[3]) # in case of 6 (1-index), cfDNA start loc will be between  delim_loc 6 (0-index) and 7 

# Create a dictionary of sites to be included in the v-plot 

flen_in_motif_tsv_fp = open(sys.argv[5])  
dna_cnt_threshold = int(sys.argv[6])
vplot_flank = int(sys.argv[7])
actual_flank = int(sys.argv[8])
min_frag_len_in_tfbs = int(sys.argv[9])
max_frag_len_in_tfbs = int(sys.argv[10])
record_tfbs_passing_th = open(sys.argv[11], "w")
print ([dna_cnt_threshold, vplot_flank, actual_flank])

to_keep_tfbs = defaultdict(lambda : False)
for line in flen_in_motif_tsv_fp: 
    # head -1 flen_count_matrix/mcf7_merged_intersect_50bp_smk_250_wMotif_er_e2_mcf7_flen_in_motif.tsv 
    # chr3	196343285	196343402	chr3:196343081-196343582|er_e2_mcf7`21620_dummy%CUTnRUN_ccaggtgccagtgaccca	60-51-79-66-35-59	6 
    # chr2	237496521	237496637	chr2:237496311-237496812|er_e2_mcf7`36418_dummy%CUTnRUN_AGGGTCATCAGAACCCA	117-44-54-38-38-57-116-150-116-97-97-117-37	13
    # get the count from the last field
    d_loc = [m.start() for m in  re.finditer("\t", line)]
    cnt = int(line[d_loc[-1]: len(line) - 1])
    if cnt >= dna_cnt_threshold:
        # select for tfbs that have at least `dna_cnt_threshold` fragments in length range
        all_lengths = [int(v) for v in line[d_loc[-2] + 1 : d_loc[-1]].split("-")] 
        to_select = False
        counter_for_th = 0
        for l in all_lengths:
            if (l >= min_frag_len_in_tfbs)  and (l<=max_frag_len_in_tfbs):
                counter_for_th = counter_for_th + 1 
        if counter_for_th >= dna_cnt_threshold:        
        # prepare the vplot compatible locus using flank
            st = int(line[d_loc[0] + 1: d_loc[1]]) - vplot_flank + actual_flank 
            en = int(line[d_loc[1] + 1: d_loc[2]]) + vplot_flank - actual_flank
            k = "\t".join([line[0 : d_loc[0]] , str(st), str(en)]) # key
            to_keep_tfbs[k] = True  
            record_tfbs_passing_th.write(line[0:d_loc[3]] + "\n")


extreme_flank_left = 100 
extreme_flank_right = -100
min_frag_len = 10000
max_frag_len = -1

vplot_matix = defaultdict(lambda : defaultdict (lambda : 0))
for line in verbose_fp:
    line_str = bytes.decode(line) 
    delim_loc = [m.start() for m in re.finditer("\t", line_str)]
    # TFBS center
    # prepare the key
    k =  line_str[0: delim_loc[2]]
    if to_keep_tfbs[k] == True:
        
        tfbs_center = int ( (int(line_str[delim_loc[0]+ 1: delim_loc[1]]) + int (line_str[delim_loc[1] + 1: delim_loc[2]]) )/2)
        cfdna_center = int ( (int(line_str[delim_loc[bed_cols]+ 1: delim_loc[bed_cols + 1]]) + int (line_str[delim_loc[bed_cols + 1] + 1: delim_loc[ bed_cols + 2]]) )/2)
        distance = cfdna_center - tfbs_center
        cfdna_frag_length = int (line_str[delim_loc[fragment_len_column - 2] + 1 : delim_loc[fragment_len_column - 1]])
        vplot_matix[cfdna_frag_length][distance] +=1
        
        if distance < extreme_flank_left: 
            extreme_flank_left = distance
        if distance > extreme_flank_right:
            extreme_flank_right = distance 
        if cfdna_frag_length < min_frag_len:
            min_frag_len = cfdna_frag_length
        if cfdna_frag_length > max_frag_len: 
            max_frag_len = cfdna_frag_length
    else:
        continue

out_fp = open(sys.argv[4], "w")
for l in range(0, max_frag_len + 1): # l mean fragment length
    # vector for l
    vec_for_l = [] 
    for f in range(extreme_flank_left, extreme_flank_right + 1): 
        vec_for_l.append(vplot_matix[l][f]) 
    to_write = "\t".join(map(str, vec_for_l))
    out_fp.write(str(l) + "\t" + to_write + "\n")

out_fp.close() 
verbose_fp.close()
record_tfbs_passing_th.close()

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

extreme_flank_left = 100 
extreme_flank_right = -100
min_frag_len = 10000
max_frag_len = -1

vplot_matix = defaultdict(lambda : defaultdict (lambda : 0))
for line in verbose_fp:
    line_str = bytes.decode(line) 
    delim_loc = [m.start() for m in re.finditer("\t", line_str)]
    # TFBS center
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

out_fp = open(sys.argv[4], "w")
for l in range(min_frag_len, max_frag_len + 1): # l mean fragment length
    # vector for l
    vec_for_l = [] 
    for f in range(extreme_flank_left, extreme_flank_right + 1): 
        vec_for_l.append(vplot_matix[l][f]) 
    to_write = "\t".join(map(str, vec_for_l))
    out_fp.write(str(l) + "\t" + to_write + "\n")

out_fp.close() 
verbose_fp.close()

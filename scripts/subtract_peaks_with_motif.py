import sys
from collections import defaultdict
import gzip
import re
import pickle

to_remove = open(sys.argv[1])
total_set = open(sys.argv[2])
out_fp = open(sys.argv[3], "w")
peak_dict = defaultdict(lambda : False)
for line in to_remove:
    # head -1 50bp_foxa1_wMotif_cnr.bed
    # chr1	18652	18766	chr1:18636-18666|FoxA1`dummy_seq`inhouse`rl_pk_%positive_ATGGTATTTACACAT	11.6122	+
    peak_key = line.split("|")[0].split("\t")[-1]
    peak_dict[peak_key] = True

# head -1 ../50bp_foxa1_e2_inhouse_cnr_peaks.bed 
# chr11	469841	469941	chr11:469876-469906|FoxA1`dummy_seq`inhouse`rl_pk_%positive	6.0	.
for line in total_set:
    peak_key = line.split("|")[0].split("\t")[-1]
    if peak_dict[peak_key] == False:
        out_fp.write(line)

out_fp.close() 
total_set.close() 
to_remove.close() 

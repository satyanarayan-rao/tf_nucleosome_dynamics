import sys
from collections import defaultdict
import re
# query file
# head -1 50bp_common_er_inh_cnr.bed
# chr1	938406	938506	chr1:938406-938506|ER`16225_dummy%CUTnRUN	1.00	-
q_fp = open(sys.argv[1])

# subject file 
# head -1 ~/softlinks/er_len_dist/flen_count_matrix/ih02_intersect_50bp_wMotifs_08_MCF7_ER_E2_peak_peaks_flen_in_motif.tsv
# chr1	11176	11276	chr1:11176-11276|ER`10891_dummy%CUTnRUN	112-140-150-67-61-59-343-163-125-109-169-173-74-150-172-82-161-59-98	19
s_fp = open(sys.argv[2])

o_fp = open(sys.argv[3], "w")

s_dict = {} 
for line in s_fp:
    d_loc = [m.start() for m in re.finditer("\t", line)] 
    k = line[0:d_loc[2]] 
    s_dict[k] = line
for line in q_fp:
    d_loc = [m.start() for m in re.finditer("\t", line)] 
    k = line[0:d_loc[2]]
    l = s_dict[k] 
    o_fp.write(l)

o_fp.close()
s_fp.close()
q_fp.close()

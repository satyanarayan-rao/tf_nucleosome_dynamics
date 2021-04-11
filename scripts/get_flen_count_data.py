import sys
from collections import defaultdict
import re

s_fp = open(sys.argv[1])
q_fp = open(sys.argv[2])
out_fp = open(sys.argv[3], "w")

s_dict = {} 
for line in s_fp:
    # example line: head -1 ../../flen_count_matrix/ih02_intersect_50bp_common_v2_flen_in_motif.tsv 
    # chr1	904715	904833	chr1:904567-904991|CTCF`74760_GCGCCCCCTGGTGGCGGAG%ENCODE	197-142-140-61-73-73-70-100	8
    d_loc = [m.start() for m in re.finditer("\t", line)] 
    k = line[0:d_loc[2]] 
    s_dict[k] = line

for line in q_fp:
    d_loc = [m.start() for m in re.finditer("\t", line)] 
    k = line[0:d_loc[2]]     
    if k in s_dict:
        to_write = s_dict[k]
        out_fp.write(to_write)
    else:
        print("{k} not found in {s}".format(k = k, s=sys.argv[1]))

out_fp.close()
s_fp.close()
q_fp.close()

import sys
from collections import defaultdict
import gzip
import re
import pickle

cfdna_map_list = sys.argv[1:len(sys.argv)] 

was_present_dict = defaultdict(lambda: False)

for f in cfdna_map_list: 
    fp = gzip.open(f)
    for line in fp: 
        # example line: zcat cfdna_sample_intersect_to_sites/ih02_on_50bp_er_e2_mcf7_th_05_1kb.bed.gz | head -1 
        # chr1	9101	11102	chr1:10051-10152|er_e2_mcf7_th_05`103609_dummy%CUTnRUN	1.0	.	chr1	9994	10123	IH02_05.pairs@1	129	+
        line_str = line.decode()
        d_loc = [m.start() for m in re.finditer("\t", line_str)] 
        cfdna_loc = line_str[d_loc[5] + 1: d_loc[8]] 
        if was_present_dict[cfdna_loc] == False:
            was_present_dict[cfdna_loc] = f
            print(line_str[d_loc[5] + 1: len(line_str)], end="") 
        elif was_present_dict[cfdna_loc] == f:
            print(line_str[d_loc[5] + 1: len(line_str)], end="") 
        elif was_present_dict[cfdna_loc] != f: 
            continue
        
    fp.close() 


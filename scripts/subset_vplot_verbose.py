import sys
from collections import defaultdict
import gzip
import re
import pickle

q_fp = open(sys.argv[1])

q_dict = defaultdict(lambda : False) 

for line in q_fp:
    d_loc = [m.start() for m in re.finditer("\t", line)]
    q_dict[line[0:d_loc[2]]] = True  
q_fp.close()
  
for line in sys.stdin:
    # example line: zless  vplot_flen_count_matrix/mcf7_merged_intersect_50bp_foxa1_e2_inhouse_cnr_peaks_total_flank_500_flen_in_motif.tsvverbos.tsv.gz | head -2
    # chr20	47602071	47603071	chr20:47602556-47602586|FoxA1`dummy_seq`inhouse`rl_pk_%positive	6	.	chr20	47602558	47602596	38	.	-
    # chr2	11502191	11503191	chr2:11502676-11502706|FoxA1`dummy_seq`inhouse`rl_pk_%positive	6	.	chr2	11503091	11503144	53	.	-
    d_loc = [m.start() for m in re.finditer("\t", line)] 
    key = line[0: d_loc[2]]
    if q_dict[key] == True:
        print(line, end = "")

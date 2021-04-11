import sys
from collections import defaultdict
import gzip
import re
import pickle

bed_dict = defaultdict (list) 
bed_in_order = [] 
bed_flag = defaultdict (lambda : False)
for line in sys.stdin:
    # example line:  zcat try_bedtools_based_mapping/pt4_on_er_th_05.bed.gz | head -2
    # chr1	845721	845822	chr1:845721-845822|er_e2_mcf7_th_05`103676_dummy%CUTnRUN	1.0	.	chr1	845763	845764	39	.	-	chr1	845744	845783	39	.	-
    # chr1	933251	933352	chr1:933251-933352|er_e2_mcf7_th_05`103680_dummy%CUTnRUN	1.0	.	chr1	933251	933252	40	.	+	chr1	933231	933271	40	.	+
    d_loc = [m.start() for m in re.finditer("\t", line)]
    key = line [0:d_loc[3]] 
    flen = line[d_loc[8] + 1: d_loc[9]]
    flen = str (int(line[d_loc[13] + 1: d_loc[14]]) - int(line[d_loc[12] + 1: d_loc[13]]))
    bed_dict[key].append(flen)

    if bed_flag[key] == False:
        bed_in_order.append(key)
        bed_flag[key] = True

for key in bed_in_order:
    # desired output: head -2 ../flen_count_matrix/pt4_merged_intersect_50bp_er_e2_mcf7_th_05_flen_in_motif.tsv
    # chr13	33581261	33581362	chr13:33581261-33581362|er_e2_mcf7_th_05`132924_dummy%CUTnRUN	42-42-37-76	4
    # chr14	31787541	31787642	chr14:31787541-31787642|er_e2_mcf7_th_05`91390_dummy%CUTnRUN	143-143-38-143-150-42	6
    
    lengths = "-".join(bed_dict[key])
    to_write = "\t".join([key, lengths, str(len(bed_dict[key])) ])
    print (to_write)

import sys
from collections import defaultdict
import gzip
import re
import pickle
q_fp = open(sys.argv[1])
s_fp = open(sys.argv[2])
out_fp = open(sys.argv[3], "w")

s_fold_dict = defaultdict(lambda : 0)
s_q_dict = defaultdict(lambda : 0)

for line in s_fp:
    #grep "905276" ~/data/encode/pu.1/bams/macs2_peaks_q_0.20/q_.20_peaks.narrowPeak 
    #chr1	905276	905493	q_.20_peak_2	95	.	7.71767	12.21043	9.55227	143
    #chr4	26905276	26905453	q_.20_peak_76652	75	.	6.47614	10.15807	7.56040	93
    # 7.71767 -> fold enrichment
    # 9.555 -> -log10(q-val)

    d_loc = [m.start() for m in re.finditer("\t", line)]
    k = line[0:d_loc[2]]
    fold = line[d_loc[5] + 1: d_loc[6]]
    q_val = line[d_loc[7] + 1 : d_loc[8]] 
    s_fold_dict[k] = fold
    s_q_dict[k] = q_val
header = "\t".join(["cl", "fold", "q_val", "sample", "peak_loc"])
out_fp.write(header + "\n")
for line in q_fp:
    # 6	chr10	100014021	100014122	chr10:100014052-100014238|pu1_peaks_q_th_point_20`14374_ggcgaagaggaagttcacct
    # 6	chr10	100076971	100077072	chr10:100076853-100077108|pu1_peaks_q_th_point_20`14378_TTTTAAGAGGAAGGAGAAAA
    d_loc = [m.start() for m in re.finditer("\t", line)] 
    cl = line[0:d_loc[0]]
    k = line[d_loc[3] + 1:].split("|")[0].replace(":", "\t").replace("-", "\t") 
    k_to_select = line[d_loc[3] + 1:].split("|")[0]
    fold = s_fold_dict[k]
    q_val = s_q_dict[k]
    to_write = "\t".join([cl, fold, q_val, sys.argv[1], k_to_select])
    out_fp.write(to_write + "\n")
    
out_fp.close() 
q_fp.close() 
s_fp.close()

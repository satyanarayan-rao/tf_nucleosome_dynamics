import sys
from collections import defaultdict
import gzip
import re
import pickle

s_fp = open(sys.argv[1])
q_fp = open(sys.argv[2]) 

s_dict = defaultdict (lambda : 0) 

for line in s_fp:
    # head -2 ../input_bed/50bp_pu1_q_30_nooverlap.bed 
    # chr1	905249	905350	chr1:905276-905493|pu1_peaks_q_th_point_30`2_gaaaagcggacgta	11.2182	+
    # chr1	906834	906935	chr1:906798-907121|pu1_peaks_q_th_point_30`3_TTAAAGGGGAAGCGACA	13.1818	-
    d_loc = [m.start() for m in re.finditer("\t", line)] 
    k = line[d_loc[2] + 1: d_loc[3]]
    s_dict[k] = line[d_loc[3] + 1: d_loc[4]] 

for line in q_fp:
    # head -2 ih02_q_30.tsv
    # 6	chr10	100014021	100014122	chr10:100014040-100014238|pu1_peaks_q_th_point_30`16343_ggcgaagaggaagttcacct
    # 6	chr10	100076971	100077072	chr10:100076853-100077108|pu1_peaks_q_th_point_30`16347_TTTTAAGAGGAAGGAGAAAA
    d_loc = [m.start() for m in re.finditer("\t", line)] 
    k = line[d_loc[3] + 1: len(line) - 1] 
    cl = line [0: d_loc[0]]
    score = s_dict[k]
    print ("{cl}\t{score}".format(cl = cl, score = score))
s_fp.close() 
q_fp.close()
   

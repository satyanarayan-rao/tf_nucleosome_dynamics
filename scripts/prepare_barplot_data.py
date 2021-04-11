import os 
import sys
import re
from collections import defaultdict
flen_cnt_fp = open(sys.argv[1])
input_bed_fp = open (sys.argv[2])
out_fp = open (sys.argv[3], "w")

cnt_dict = defaultdict(lambda : 0) 

for line in flen_cnt_fp: 
    d_loc = [m.start() for m  in re.finditer("\t", line)]
    cnt_dict[line[d_loc[2]+1:d_loc[3]]] = int (line[d_loc[-1]:len(line)])
for line in input_bed_fp: 
    d_loc = [m.start() for m  in re.finditer("\t", line)]
    k = line[d_loc[2]+1:d_loc[3]]
    to_write = "\t".join([k, str(cnt_dict[k])])
    out_fp.write (to_write + "\n")
out_fp.close() 
input_bed_fp.close()
flen_cnt_fp.close() 

import os
import sys
import re
from collections import defaultdict
import gzip
inp_fp = open(sys.argv[1])
cutoff = float(sys.argv[2])
out_fp = open(sys.argv[3], "w")
for line in inp_fp:
    line_length = len(line)
    d_loc = [m.start() for m in re.finditer("\t", line)]
    peak_distance = abs(float(line[d_loc[-1] + 1: len(line) - 1]))
    if peak_distance > cutoff: 
        to_write = line[0: d_loc[3]] + "`" + ">" + str(cutoff) + line[d_loc[3]: line_length] 
        out_fp.write(to_write) 
    elif  (peak_distance > 0) and  (peak_distance<=100): 
       to_write = line[0: d_loc[3]] + "`" + "0-100" + line[d_loc[3]: line_length] 
       out_fp.write(to_write)
    elif  (peak_distance > 100) and  (peak_distance<=200): 
       to_write = line[0: d_loc[3]] + "`" + "101-200" + line[d_loc[3]: line_length] 
       out_fp.write(to_write)
    elif  (peak_distance > 200) and (peak_distance<=400): 
       to_write = line[0: d_loc[3]] + "`" + "201-400" + line[d_loc[3]: line_length] 
       out_fp.write(to_write)
    elif  (peak_distance > 400) and  (peak_distance<=800): 
       to_write = line[0: d_loc[3]] + "`" + "401-800" + line[d_loc[3]: line_length] 
       out_fp.write(to_write)
    elif  (peak_distance > 801) and  (peak_distance<=cutoff): 
        to_write = line[0: d_loc[3]] + "`" + "801-" + str(cutoff)  + line[d_loc[3]: line_length] 
        out_fp.write(to_write)
out_fp.close() 
inp_fp.close()

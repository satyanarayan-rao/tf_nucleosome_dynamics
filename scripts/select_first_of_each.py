import sys
import re

inp_fp = open(sys.argv[1])

prev = "somethig"

for line in inp_fp:
    if len(line) <5: 
        continue
    else:
        d_loc = [m.start() for m in re.finditer("\t", line)] 
        peak_id = line[d_loc[1] + 1: d_loc[2]] # mind that there are two tabs  
        if peak_id !=prev: 
            print(line, end='') 
            prev=peak_id
        else:
            continue

import os
import sys
from collections import defaultdict
import re
input_fp = open(sys.argv[1])
total_dict = defaultdict(lambda : 0 ) 
out_fp = open("segway_genome_percent.tsv", "w")
g_size=3000000000
for line in input_fp:
    delim_loc = [m.start() for m in re.finditer("\t", line)]
    start=int(line[delim_loc[0]+1: delim_loc[1]])
    end=int(line[delim_loc[1]+1: delim_loc[2]])
    label = line[delim_loc[2] +1: delim_loc[3]]
    total_dict[label] += end - start

for label in total_dict: 
    val = total_dict[label]
    per = round(val/g_size, 5)
    to_write = "\t".join([label, str(val), str(per)])
    out_fp.write(to_write + "\n")

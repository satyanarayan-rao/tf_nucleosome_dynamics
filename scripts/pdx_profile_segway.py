import os
import sys
import gzip
import re
import pickle
from collections import defaultdict

input_fp = open(sys.argv[1])
out_fp = open ("chrom_stats.tsv", "w")

total_mapped_reads = 0
read_dict = defaultdict(list)

for line in input_fp: 
    delim_loc = [m.start() for m in re.finditer("\t", line)]
    kk = line[delim_loc[8] + 1: delim_loc[11]]
    read_dict[kk].append(line[delim_loc[2]+1:delim_loc[3]])
    total_mapped_reads +=1

label_dict = defaultdict(lambda: 0 )
for read in read_dict: 
    for label in read_dict[read]: 
        label_dict[label] +=1

for label in label_dict:
    val = label_dict[label]
    per = round(val/total_mapped_reads, 2)
    to_write = "\t".join([label, str(val), str(per)])
    out_fp.write(to_write + "\n")

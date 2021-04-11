import os 
import sys
import pickle
import re
genome_cov_fp = open(sys.argv[1])
to_normalize = open(sys.argv[2])
out_fp = open(sys.argv[3], "w") 

genome_bp = 0
genome_per = 0 
for line in genome_cov_fp: 
    line_items = line.split()
    if "Quiescent" in line: 
        genome_bp += int(line_items[1]) 
        genome_per += float(line_items[1])
    else:
        out_fp.write(line)

out_fp.write("\t".join(["Quiescent", str(genome_bp), str(genome_per)]) + "\n")

import os
import sys 
import pickle
import re
# argv[1]: neg ctl bed
# argv[2]: chip motif pkl
# argv[3]: output file

input_fp = open (sys.argv[1])
chip_motif_dict = pickle.load(open(sys.argv[2], "rb"))
kk = list(chip_motif_dict.keys()) 
output_fp = open(sys.argv[3], "w")
for line in input_fp:
    delim_loc = [m.start() for m in re.finditer("\t", line)]
    key_to_find = line[0:delim_loc[2]]
    if key_to_find not in chip_motif_dict:
        output_fp.write(line)
    else:
        continue

import os
import sys
import pickle
import re
# argv[1]: chip bed file
# argv[2]: chip bed pkl file

input_fp = open(sys.argv[1])
output_fp = open(sys.argv[2], "wb")
motif_dict = {} 
for line in input_fp: 
    delim_loc =  [m.start() for m in re.finditer("\t", line)]
    key_str = line[0:delim_loc[2]]
    motif_dict[key_str] = True
    
pickle.dump (motif_dict, output_fp)
output_fp.close() 
input_fp.close()

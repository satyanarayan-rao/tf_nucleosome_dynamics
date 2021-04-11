import os
import sys 
import re
import pickle
input_fp = open(sys.argv[1])
cell_line_dict = pickle.load(open(sys.argv[2], "rb"))
out_fp = open (sys.argv[3], "w")

for line in input_fp:
    delim_loc = [m.start() for m in re.finditer("\t", line)]
    key = line[0:delim_loc[0]] + ":" + \
          line[delim_loc[0]+1:delim_loc[1]] + "-" + \
          line[delim_loc[1]+1:delim_loc[2]]
    cell_lines = "&".join(cell_line_dict[key])
    out_line = line[0:delim_loc[2] + 1]
    appended_name = line[delim_loc[2] + 1: delim_loc[3]] + "`" + cell_lines
    out_line = out_line + appended_name + line[delim_loc[3]: len(line)]
    out_fp.write(out_line)

out_fp.close()
           

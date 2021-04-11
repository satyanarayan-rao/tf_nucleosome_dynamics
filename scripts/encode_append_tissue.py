import os 
import sys
import pickle
import re 
input_fp = open (sys.argv[1])
tissue_cell_line_string = sys.argv[2]
out_fp = open(sys.argv[3], "w")

tissue_list = tissue_cell_line_string.split("^")
cell_line_dict = {}
for i in tissue_list: 
    cell_line = i.split("@")[0]
    cell_line_dict[cell_line] = i.split("@")[1]
for line in input_fp:
    delim_loc = [m.start() for m in re.finditer("\t", line)]
    name = line[delim_loc[2]+1:delim_loc[3]] 
    cell_lines = name.split("`")[2].split("&")
    tmp_tissue_list = []
    for c in cell_lines: 
        tmp_tissue_list.append(cell_line_dict[c])
    set_tmp_tissue_list = set(tmp_tissue_list) 
    if len(set_tmp_tissue_list) == 1: 
        out_line = line[0:delim_loc[2]+1] 
        name_field = line[delim_loc[2]+1:delim_loc[3]] + "`" + tmp_tissue_list[0]
        out_line = out_line + name_field + line[delim_loc[3]:len(line)]
        out_fp.write(out_line)
    else:
        out_line = line[0:delim_loc[2]+1] 
        name_field = line[delim_loc[2]+1:delim_loc[3]] + "`" + "common"
        out_line = out_line + name_field + line[delim_loc[3]:len(line)]
        out_fp.write(out_line)
        
out_fp.close() 


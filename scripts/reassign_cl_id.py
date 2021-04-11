import sys
import os
import pandas as pd
tsv_fp = open(sys.argv[1])
reassign_cl_id_fp = pd.read_csv(sys.argv[2],
                     sep="\s+", comment = "#", header = "infer")
out_fp = open(sys.argv[3], "w")
actual_list = list(reassign_cl_id_fp.actual)
assigned_list = list(reassign_cl_id_fp.assigned) 
cl_reassign_dict = {} 
for actual, assigned in zip(actual_list, assigned_list):
    cl_reassign_dict[str(actual)] = str(assigned)

for line in tsv_fp:
    line_items = line.strip().split() 
    line_items[0] = cl_reassign_dict[line_items[0]] 
    to_write = "\t".join(line_items) 
    out_fp.write(to_write + "\n")

out_fp.close() 
tsv_fp.close()

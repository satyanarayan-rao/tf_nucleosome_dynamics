import os
import sys
import pickle
import re
from collections import defaultdict
encode_fp = open(sys.argv[1])
carroll_fp = open(sys.argv[2])
tmp_encode_file = "tmp/" + sys.argv[1].split("/")[-1]
tmp_carroll_file = "tmp/" + sys.argv[2].split("/")[-1]

tmp_encode_fp = open(tmp_encode_file, "w")
tmp_carroll_fp = open(tmp_carroll_file, "w")

motif_dict = defaultdict(list)
for line in encode_fp:
    delim_loc = [m.start() for m in re.finditer("\t", line)] 
    key = line[0:delim_loc[2]]
    motif_dict[key].append("ENCODE") 

for line in carroll_fp: 
    delim_loc = [m.start() for m in re.finditer("\t", line)] 
    key = line[0:delim_loc[2]]
    motif_dict[key].append("Carroll") 

encode_fp.close() 
carroll_fp.close()

encode_fp = open(sys.argv[1])
carroll_fp = open(sys.argv[2])

for line in encode_fp:
    delim_loc = [m.start() for m in re.finditer("\t", line)] 
    key = line[0:delim_loc[2]]
    test_for_common = motif_dict[key]
    set_test_for_common  = set(test_for_common)
    if len (set_test_for_common) == 1: # unique to either ENCODE or Carroll
        to_write = line[0: delim_loc[2] + 1] + line[delim_loc[2] + 1: delim_loc[3]] + "`Unique_ENCODE" + line[delim_loc[3]: len(line)]
        tmp_encode_fp.write(to_write)
    else:
        to_write = line[0: delim_loc[2] + 1] + line[delim_loc[2] + 1: delim_loc[3]] + "`Common_in_ENCODE_and_Carroll" + line[delim_loc[3]: len(line)]
        tmp_encode_fp.write(to_write)
tmp_encode_fp.close()

for line in carroll_fp:
    delim_loc = [m.start() for m in re.finditer("\t", line)] 
    key = line[0:delim_loc[2]]
    test_for_common = motif_dict[key]
    set_test_for_common  = set(test_for_common)
    if len (set_test_for_common) == 1: # unique to either ENCODE or Carroll
        to_write = line[0: delim_loc[2] + 1] + line[delim_loc[2] + 1: delim_loc[3]] + "`Unique_Carroll" + line[delim_loc[3]: len(line)]
        tmp_carroll_fp.write(to_write)
    else:
        to_write = line[0: delim_loc[2] + 1] + line[delim_loc[2] + 1: delim_loc[3]] + "`Common_in_ENCODE_and_Carroll" + line[delim_loc[3]: len(line)]
        tmp_carroll_fp.write(to_write)
tmp_carroll_fp.close()

# cat both tmp files and do sort to select for in

cmd_to_run = "cat {f1} {f2} | sort -k1,1 -k2,2n -k3,3n --stable --unique > {out_file}".format(f1=tmp_encode_file, f2 = tmp_carroll_file, out_file = sys.argv[3])
print(cmd_to_run)
os.system (cmd_to_run)

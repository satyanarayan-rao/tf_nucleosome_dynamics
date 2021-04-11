import os
import sys
import re
from collections import defaultdict

# from_ih02_130_180.bed
# mcf7.bed
# from  /beevol/home/satyanarr/softlinks/ctcf_len_dist/heatmap_one_system_to_other
query_fp = open(sys.argv[1])
subject_fp = open(sys.argv[2])

# example line qery_fp
# chr10	100239910	100242028	chr10:100240692-100241087|peak8269_CTGTCTCTAGGGGGAGCAC%Carroll^6
# example line subject_fp 
# chr10	100018957	100021075	chr10:100019755-100020234|CTCF`159918_GAGACTGCAGGGGGAGGCA%ENCODE^6

subject_dict = defaultdict(lambda : "NA")
count_dict = defaultdict(lambda : 0)
for line in subject_fp:
    d_loc = [m.start() for m in re.finditer("\t", line)]
    cl_id = line[d_loc[-1]: len(line) - 1].split("^")[-1]  
    subject_dict[line[0:d_loc[2]]] = cl_id
    count_dict[cl_id] +=1

reflection_dict = defaultdict(lambda : defaultdict (lambda : 0))
for line in query_fp:
    d_loc = [m.start() for m in re.finditer("\t", line)]
    cl_id = line[d_loc[-1]: len(line) - 1].split("^")[-1]
    sub_cl = subject_dict[line[0:d_loc[2]]]
    if sub_cl != "NA":
        reflection_dict[cl_id][sub_cl] +=1

    #print (line[0:d_loc[2]] + "\t" + cl_id + "\t" + subject_dict[line[0:d_loc[2]]])
order = list (map(str, range(1,7)))    
last_line = [] 
for k1 in order:
    print (k1, end = "\t")
    total_count = 0 
    for k2 in order:
        print (str(reflection_dict[k1][k2]), end = "\t")
        total_count += reflection_dict[k2][k1]
    last_line.append(count_dict[k1] - total_count)
    print ()
print ("\t".join([str(len(order) + 1), "\t".join(map(str, last_line))]))

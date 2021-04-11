import sys
from collections import defaultdict
import re
from itertools import zip_longest
bed_fp = open(sys.argv[1]) 
cnt_th = int(sys.argv[2])
out_fp = open(sys.argv[3], "w")
list_of_files = sys.argv[4:]
total_files = len(list_of_files)


bed_dict_fourth_field = {} 
bed_dict_first_three = {} 
out_cnt_fnames = [l.split("/")[-1].split("_intersect_")[0] for l in list_of_files] 
print (out_cnt_fnames)
out_cnt_fp_list = [open("common_" + out_cnt_fnames[i] + ".tsv", "w") for i in range(total_files)]

fp_list = [open(f) for f in list_of_files]
for line in bed_fp:
     d_loc = [m.start() for m in re.finditer("\t", line)] 
     key_first_three = line[0:d_loc[2]] 
     #key_four = line[d_loc[2]+1: d_loc[3]] 
     key_four = line[0: d_loc[3]] 
     bed_dict_fourth_field[key_four] = line
     bed_dict_first_three[key_first_three] = line

cnt = 0 
list_of_dict = [defaultdict(lambda : 0) for i in range(len(list_of_files))]

list_of_cfdna_cnt_dict = [{} for i in range(total_files) ] 
for f in fp_list:
    
    for line in f:
        d_loc = [m.start() for m in re.finditer("\t", line)]
        #k = line[d_loc[2] +1: d_loc[3]]
        k = line[0: d_loc[3]]
        v = int(line[d_loc[-1]+1:]) 
        list_of_dict[cnt][k] = v
        list_of_cfdna_cnt_dict[cnt][k] = line 
        
    cnt +=1     
    f.close()

for k in list_of_dict[0]:
    is_gt_th = True
    for j in range(total_files):
        if list_of_dict[j][k] < cnt_th:
            is_gt_th = False
    if is_gt_th == True:
        out_fp.write(bed_dict_fourth_field[k])
        for q in range(total_files):
             out_cnt_fp_list[q].write(list_of_cfdna_cnt_dict[q][k])

out_fp.close()

for q in range(total_files):
    out_cnt_fp_list[q].close()








#cnt = 0 
#for lines in zip_longest(*fp_list, fillvalue=''):
#    # example line: head -1 flen_count_matrix/ih02_intersect_50bp_excluded_flen_in_motif.tsv 
#    # chr1	154193917	154194017	chr1:154193917-154194017|ER`2_dummy%CUTnRUN	177-243-139-132-132-174-177-157-166-161-161-177-162-163-95-63-74-71-71	19
#    #print(lines)
#    tab_of_all = [[m.start() for m in re.finditer("\t", line)] for line in lines]
#    #print(tab_of_all[0][-1])
#    #print (tab_of_all)
#    print(lines[1], end = '') 
#    read_cnt_in_all = [int(lines[i][tab_of_all[i][-1]: len(lines[i]) - 1]) for i in range(total_files)]   
#    key_in_all = [lines[i][tab_of_all[i][2]: tab_of_all[i][3]] for i in range(total_files)]
#    for i, k, v in  zip(range(total_files), read_cnt_in_all, key_in_all): 
#        list_of_dict[i][k] = v   
#
#
#for fp in fp_list:
#    fp.close()

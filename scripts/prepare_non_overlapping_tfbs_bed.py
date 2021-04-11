import os
import sys
from collections import defaultdict
import re

inp_fp = open(sys.argv[1]) 
th_flank = int(sys.argv[2])
out_fp = open(sys.argv[3], "w")


locus_to_line_map = defaultdict(lambda : False)
locus_to_score_map = defaultdict(lambda : False)
overlapped_and_best = defaultdict(lambda : False)
discard_lines = defaultdict(lambda : False)
score_boundary_st = defaultdict(lambda : False)
score_boundary_en = defaultdict(lambda : False)
line_dict = {} 
cnt = 1 

to_print_dict = {} # defaultdict (lambda : False)
for line in inp_fp:
    line_dict[cnt] = line
    cnt +=1

inp_fp.close()

inp_fp = open(sys.argv[1]) 
cnt = 1 
for line in inp_fp:
    # example lines:  head -5 input_bed/50bp_common_v2.bed
    #chr1	904715	904833	chr1:904567-904991|CTCF`74760_GCGCCCCCTGGTGGCGGAG%ENCODE	9.0	-
    #chr1	904717	904835	chr1:904567-904991|CTCF`74760_CCGCCACCAGGGGGCGCCA%ENCODE	25.623	+
    #chr1	908139	908257	chr1:908109-908339|CTCF`74761_TGTCCACGAGGAGGACCGG%ENCODE	7.93443	-
    #chr1	908177	908295	chr1:908109-908339|CTCF`74761_ggcccacaggagggcgggc%ENCODE	10.5246	+
    #chr1	912704	912822	chr1:912712-913240|CTCF`74762_TTGGCAGAAGGTGGCTCTG%ENCODE	11.540999999999999	+
    # Condition : bed file should have unique coordinates (1-3 should be unique)
    d_loc = [m.start() for m in re.finditer("\t", line)] 
    chrom = line[0:d_loc[0]] 
    st = int(line[d_loc[0] + 1: d_loc[1] ])
    en = int(line[d_loc[1] + 1: d_loc[2]])
    key = line[0:d_loc[2]] 
    score = float(line[d_loc[3]+ 1: d_loc[4]]) 
    #print([key,score])
    #print ([st - th_flank, en + th_flank])
    for bp in range(st - th_flank, en + th_flank):
        bp_key = chrom + "-" + str(bp)
        if locus_to_line_map[bp_key] == False:
            locus_to_line_map[bp_key] = cnt 
            locus_to_score_map[bp_key] = score
            to_print_dict[bp_key] = cnt
            score_boundary_st[bp_key]  = st - th_flank
            score_boundary_en[bp_key]  = en + th_flank
            #print([cnt, bp_key, locus_to_score_map[bp_key]])
        else:
            #print ([bp_key, score])
            if locus_to_score_map[bp_key] < score: 
                  #print ([bp_key, to_print_dict[bp_key], locus_to_score_map[bp_key], cnt, score])
                  locus_to_score_map[bp_key] = score
                  #print (to_print_dict[bp_key])
                  discard_lines[to_print_dict[bp_key]] = True
                  to_print_dict[bp_key] = cnt 
                  #print (cnt)
                  #sys.exit(-1)
                  overlapped_and_best[cnt] = True # bigger score came first, so discarding the second
                  max_en = max(en + th_flank, score_boundary_en[bp_key])
                  min_st = min(st - th_flank, score_boundary_st[bp_key])
                  #print ([min_st, max_en])
                  for pp in range(min_st, max_en): 
                      pk = chrom + "-" + str(pp)
                      locus_to_score_map[pk] = score
                      to_print_dict[pk] = cnt # over-write the line index with new line number 
                      #print ([pk, to_print_dict[pk]])
                  #print(["#####\t", to_print_dict[bp_key]])
                  #print (["@@@\t", bp_key, to_print_dict[bp_key], locus_to_score_map[bp_key], cnt, score])
                  break
            else:
                discard_lines[cnt] = True
                break
             
    cnt +=1         

already_printed_dict = defaultdict(lambda : False)

for k in range(1, cnt): 
    if discard_lines[k] == False:
        out_fp.write(line_dict[k])
    else:
        print(line_dict[k].strip()) 
        continue

out_fp.close() 

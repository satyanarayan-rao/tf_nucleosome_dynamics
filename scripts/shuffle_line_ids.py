import os
import sys
import random
SEED = 731 
tot_lines = int (sys.argv[1])
output_fp = open(sys.argv[2], "w")
line_list = list(range(1, tot_lines + 1 ))
random.seed(SEED)
random.shuffle(line_list)
for l in line_list: 
    output_fp.write(str(l) + "\n")

output_fp.close() 


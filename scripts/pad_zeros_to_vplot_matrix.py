import sys
from collections import defaultdict
import re

inp_fp = open (sys.argv[1]) 
out_fp = open (sys.argv[2], "w")

first_line = inp_fp.readline().split()
min_frag = int(first_line[0])
xra = len(first_line) - 1
for i in range(min_frag):
    zeros = str(i) + "\t" + "\t".join(["0"]*xra) 
    out_fp.write(zeros + "\n")

inp_fp.seek(0)
for line in inp_fp:
    out_fp.write(line)
out_fp.close()
inp_fp.close()  
   

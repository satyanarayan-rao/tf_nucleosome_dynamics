import sys 
from collections import defaultdict 
import re 
import gzip
input_fp = gzip.open (sys.argv[1])
out_fp1 = open (sys.argv[2], "w")
out_fp2 = open (sys.argv[3], "w")
len_dict = defaultdict (lambda : 0 )
cnt = 0 
pattern = re.compile(bytes('\t', encoding = 'ascii'))
for line in input_fp:
    all_delims = [m.start() for m in re.finditer(pattern, line)] 
    st = int(line[all_delims[0]:all_delims[1]])
    end = int(line[all_delims[1]: all_delims[2]]) 
    length = end - st
    len_dict [length]  +=1 
    cnt +=1
 
all_keys = sorted (list (len_dict.keys()))
for k in all_keys: 
    to_write = str (k) + "\t" + str (round (len_dict[k]/cnt, 4))
    out_fp1.write (to_write + "\n") 
    to_write = str (k) + "\t" + str (len_dict[k])
    out_fp2.write (to_write + "\n") 

input_fp.close() 
out_fp1.close() 
out_fp2.close() 

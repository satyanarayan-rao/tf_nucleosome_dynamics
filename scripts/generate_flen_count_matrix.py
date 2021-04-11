import os
import sys
import pickle
from collections import defaultdict
import gzip
import re

def initiate(max_length):
    dd = {} 
    for i in range (max_length + 1): 
        dd[i] = 0 
    return dd
"""
chr1    1421146 1421560 chr1:1421240-1421475|FOXA1`778251_ACGCTGTTTACTCAG%ENCODE        15.4082 +       chr1    1421111 1421278 167    .       +
"""
input_fp = gzip.open(sys.argv[1])
out_fp = open (sys.argv[2], "w")
out_motif_fp = open(sys.argv[3], "w")

motif_flen_counts = defaultdict(lambda : [0]*701)
flen_in_motifs = defaultdict(list)
pattern = re.compile(bytes("\t", encoding = 'ascii'))

for line in input_fp: 
    all_delim_start_loci = [m.start() for m in re.finditer(pattern, line)] 
    # since the bed file contains no repeating motifs, I am going to use first three fields as the key
    tfbs_motif = bytes.decode(line[0:all_delim_start_loci[3]])
    flen = int(bytes.decode(line[all_delim_start_loci[8]+1: all_delim_start_loci[9]]))
    #print([tfbs_motif, flen])
    motif_flen_counts[tfbs_motif][flen] += 1
    flen_in_motifs[tfbs_motif].append(flen)
    

for k in motif_flen_counts.keys(): 
    count_vec_str = "\t".join (map(str, motif_flen_counts[k]))
    to_write = "\t".join([k, count_vec_str])
    out_fp.write(to_write + "\n")
    flen_in_motif_str = "-".join(map(str, flen_in_motifs[k]))
    to_write1 = "\t".join([k, flen_in_motif_str, str (len(flen_in_motifs[k]))])
    out_motif_fp.write(to_write1 + "\n")

out_motif_fp.close()
out_fp.close()
input_fp.close()

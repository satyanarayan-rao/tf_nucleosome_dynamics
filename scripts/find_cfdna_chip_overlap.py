import os
import sys
from collections import defaultdict
# sys.argv[1]: chip_kmeans_row_order
# sys.argv[2]: cfdna clsuter
# sys.argv[3]: overlap output file
chip_kmean_fp = open(sys.argv[1])
cfdna_fp = open(sys.argv[2])
overlap_fp = open (sys.argv[3], "w")
overlap_fp.write("chr_loc\tchip_kmean_cl_id\tcfdna_cl_id\n")
# first creat cfdna dictionary
cfdna_dict = defaultdict(lambda : "not.found")
for line in cfdna_fp:
    line_items = line.strip().split()
    adjusted_chr_start = str(int(line_items[2]) - 250)
    adjusted_chr_end = str(int(line_items[3]) + 250)
    key = "@".join([line_items[1], adjusted_chr_start, 
                    adjusted_chr_end, line_items[4]])
    cfdna_dict[key] = line_items[0]
    #print (key)
    #break

for line in chip_kmean_fp:
    line_items = line.strip().split() 
    key_to_find = line_items[0].split("^")[0]
    #print (key_to_find) 
    #break
    to_write = "\t".join([key_to_find, line_items[1],
                          cfdna_dict[key_to_find]])
    overlap_fp.write(to_write + "\n")

overlap_fp.close()
cfdna_fp.close()
chip_kmean_fp.close()



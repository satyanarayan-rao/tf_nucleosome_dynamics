import sys
from collections import defaultdict
import re
## Subtracting `b` from `a`
## for example, to define mcf7_specific, a = mcf7 reassigned bed, b = ih02 reassigned bed
a_file_list = sys.argv[1].split()
a_cl_list = sys.argv[2].split() 
b_file_list = sys.argv[3].split()
b_cl_list = sys.argv[4].split()


a_sites = defaultdict (lambda: False)
b_sites = defaultdict (lambda: False)

sig_cl_dict = defaultdict(lambda: defaultdict (lambda : False))

# For both a and b create a dictionary of significant clusters
for cl_list, fname in zip(a_cl_list, a_file_list):
    cl_ids = cl_list.split("@")
    for cl in cl_ids:
        sig_cl_dict[fname][cl] = True
for cl_list, fname in zip(b_cl_list, b_file_list):
    cl_ids = cl_list.split("@")
    for cl in cl_ids:
        sig_cl_dict[fname][cl] = True
### 

### Create a dict of significant sites in `a`
a_sig_sites = defaultdict(lambda : False) 
a_sig_bed_lines = defaultdict()
for fname in a_file_list: 
    fp = open(fname)
    for line in fp:
        d_loc = [m.start() for m in re.finditer("\t", line)]
        cl_id = line.strip().split("^")[-1]
        site_loc = line[0:d_loc[2]] 
        if sig_cl_dict[fname][cl_id] == True: # check if cl_id is in significant category
            a_sig_sites[site_loc] = True  
         
        
    fp.close() 

b_sig_sites = defaultdict(lambda : False) 
b_sig_bed_lines = defaultdict()
for fname in b_file_list: 
    fp = open(fname)
    for line in fp:
        d_loc = [m.start() for m in re.finditer("\t", line)]
        cl_id = line.strip().split("^")[-1]
        site_loc = line[0:d_loc[2]] 
        if sig_cl_dict[fname][cl_id] == True: # check if cl_id is in significant category
            b_sig_sites[site_loc] = True
         
        
    fp.close()

## identify common and unique in a by subtracting b
specific_bed_fp = open(sys.argv[5], "w")
common_fp = open(sys.argv[6], "w") 
for site in a_sig_sites:
    if b_sig_sites[site] == True:
        common_fp.write(site + "\n") 
    else:
        specific_bed_fp.write(site + "\n")

common_fp.close() 
specific_bed_fp.close()

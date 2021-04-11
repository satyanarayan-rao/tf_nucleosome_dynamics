import os 
import sys
import re
import numpy as np
# sys.argv[1] = input tsv ; check tcga_atac_seq_analysis/tfbs_chip_and_mcf7_neg_ctl_mapped_to_er_pos_her2_neg_high_prolif.tsv
# sys.argv[2] = output file for boxplot

input_fp = open(sys.argv[1])
output_fp = open(sys.argv[2], "w")
for line in input_fp:
    delim_loc = [m.start() for m in re.finditer("\t", line)]
    tfbs_class = line[delim_loc[2] + 1: delim_loc[3]].split("%")[-1]
    atac_score = line[delim_loc[9] + 1: delim_loc[10]]
    raw_values = list(map(int, line[delim_loc[10] + 1:].split()))
    mean_of_patient_values = np.mean(raw_values)
    tfbs = line[delim_loc[2] + 1: delim_loc[3]]
    # first write single value along with tfbs
    to_write_atac_score = tfbs + "\t" + tfbs_class + "\t" + \
                          atac_score + "\t" + "atac_score" + "\t" + \
                          "atac_score-" + tfbs_class 
    output_fp.write(to_write_atac_score + "\n")
    
    to_write_mean_score = tfbs + "\t" + tfbs_class + "\t" + \
                          str(mean_of_patient_values) + "\t" + "mean_atac_score" + "\t" + \
                          "mean_atac_score-" + tfbs_class 
    output_fp.write(to_write_mean_score + "\n")
    for v in raw_values: 
        to_write_raw_atac = tfbs + "\t" + tfbs_class + "\t" + \
                              str(v) + "\t" + "raw_atac_score" + "\t" + \
                              "raw_atac_score-" + tfbs_class 
        output_fp.write(to_write_raw_atac + "\n") 

output_fp.close()

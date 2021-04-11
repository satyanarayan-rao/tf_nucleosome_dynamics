import sys
from collections import defaultdict
import gzip
import re
import pickle
# zcat vplot_flen_count_matrix/ih02_intersect_50bp_unique_carroll_encode_chrY_removed_non_overlapping_total_flank_500_flen_in_motif.tsvverbos.tsv.gz | head -1
# chr1	90919	91937	chr1:91156-91580|CTCF`74745_GTGGCACCAGGTGGCAGCA%ENCODE	16.2951	+	chr1	90838	91006	IH02_00.pairs@15152	168	+
inp_fp = gzip.open(sys.argv[1])

# head -1 fragment_bw_map_manual/ih02_intersect_50bp_unique_carroll_encode_chrY_removed_non_overlapping_dmatrix_flen_min_35_max_250_least_read_5_bw_3_nclust_6_max_iter_250_slop_1000_fsize_0_79_chip_GM_ENCFF578TBN_ref_cluster_6_reassigned_bed.tsv 
# 6	chr10	100240910	100241028	chr10:100240692-100241087|peak8269_CTGTCTCTAGGGGGAGCAC%Carroll
bed_fp = open(sys.argv[2])
out_fp = gzip.open(sys.argv[3], "wb")
flank = int(sys.argv[4]) # 500
actual_flank = int(sys.argv[5]) # 50 bp
bed_id_to_cl_map = defaultdict(lambda : "NA")

# create a map of TFBS +-50 to cl id
for line in bed_fp:
    d_loc = [m.start() for m in re.finditer("\t", line)]
    cl_id = line[0:d_loc[0]]
    bed_id = line[d_loc[0] + 1: d_loc[3]] 
    bed_id_to_cl_map[bed_id] = cl_id 
for line in inp_fp:
    line_str = bytes.decode(line) 
    d_loc = [m.start() for m in re.finditer("\t", line_str)]
    chrom = line_str[0: d_loc[0]]
    st = int(line_str[d_loc[0] + 1: d_loc[1]]) + flank - actual_flank
    en = int(line_str[d_loc[1] + 1: d_loc[2]])  - flank + actual_flank
    k = "\t".join ([chrom, str(st), str(en)]) 
    cl_id = bed_id_to_cl_map[k]
    to_write = bytes(cl_id +"\t", encoding="ascii") + line
    out_fp.write(to_write)

out_fp.close()
inp_fp.close() 
bed_fp.close()

import sys
import re

inp_fp = open(sys.argv[1])
out_fp = open(sys.argv[2], "w")
slop_bp = sys.argv[3]

# chr1	18490	18606	chr1:18461-18862|er_e2_mcf7`41659_dummy%CUTnRUN_AAGGATAATCTGACCTG	-0.46969700000000003	+
# chr1	186009	186125	chr1:185941-186342|er_e2_mcf7`41665_dummy%CUTnRUN_CCTGTCAGTTTGAGCTG	-1.71212	+
# chr1	189013	189129	chr1:188881-189282|er_e2_mcf7`41666_dummy%CUTnRUN_AAGGATAATCTGACCTG	-0.46969700000000003	+
# chr1	199740	199856	chr1:199691-200092|er_e2_mcf7`41668_dummy%CUTnRUN_AGGGTCGTGGGGGCCTG	2.16667	-
# chr1	938129	938245	chr1:938091-938492|er_e2_mcf7`41698_dummy%CUTnRUN_GCGGTCGGCAGGACCCA	0.712121	+
# chr1	1073765	1073881	chr1:1073731-1074132|er_e2_mcf7`41700_dummy%CUTnRUN_GTGGTCACAGTGACCTG	11.8485	-
# chr1	1079529	1079645	chr1:1079461-1079862|er_e2_mcf7`41702_dummy%CUTnRUN_GAGGTCACTCTGGGCTC	4.65152	+
# chr1	1080193	1080309	chr1:1080121-1080522|er_e2_mcf7`41703_dummy%CUTnRUN_AGGGACACACTGAGCTG	1.98485	-
# chr1	1080336	1080452	chr1:1080121-1080522|er_e2_mcf7`41703_dummy%CUTnRUN_GGTGTCAGGATGACCCT	16.8333	+
# chr1	1250123	1250239	chr1:1250151-1250552|er_e2_mcf7`41706_dummy%CUTnRUN_tgggtgggcctgaccct	1.0303	+

out_fp.write("\t".join(["peak_to_motif_distance", "slop", "strand"]) + "\n")

for line in inp_fp:
    d_loc = [m.start() for m in re.finditer("\t", line)]
    name_field = line[d_loc[2] + 1:d_loc[3]] 
    chr_loc = re.split("-|:", name_field[0:name_field.find("|")]) 
    peak_center = int ( (int (chr_loc[1]) + int (chr_loc[2]))/2)  
    motif_center = int ( (int (line[d_loc[0] + 1: d_loc[1]]) + int (line[d_loc[1] + 1: d_loc[2]]) )/2) 
    dist = motif_center - peak_center 

    strand = line[d_loc[4] + 1: len(line) - 1]
    to_write = "\t".join([str(dist), "flank_cnr_"+slop_bp, strand])
    out_fp.write(to_write + "\n")

out_fp.close()
inp_fp.close()

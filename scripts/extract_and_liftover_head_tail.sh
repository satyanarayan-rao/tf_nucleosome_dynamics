#!/bin/bash
# $1: chiapet intersected with cnr wmotif peaks
# $2: liftover chain file
# $3: output file 

# head -2 er_chia_pet_data/intersect_with_cnr/ihh015f_and_ihm001f_on_er_wmotif.bed
# chr1	1073396	1073515	chr1:1073161-1073662|er_e2_mcf7`41699_dummy%CUTnRUN_TGGTCAGGGCACCCTCAGCC	12.9113	+	chr1	1071187	1075654	chr1:1004671-1006376|head_tail
# chr1	1073396	1073515	chr1:1073161-1073662|er_e2_mcf7`41699_dummy%CUTnRUN_TGGTCAGGGCACCCTCAGCC	12.9113	+	chr1	1071928	1074971	chr1:1004805-1006952|head_tail

awk -F'[\t|_]' '{print $(NF-2)"\t"$NF}' $1 | awk -F '[-:]' '{OFS="\t"} {print $1, $2, $3, $4}' > ${3}.tmp

liftOver ${3}.tmp $2 ${3} ${3}.tmp2

# cleanup
rm ${3}.tmp ${3}.tmp2


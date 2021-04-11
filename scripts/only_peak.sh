awk -F'|' '{print $1}' $1  | awk '{print $NF}' | awk -F'[-:]' '{print $1"\t"$2"\t"$3}' | sort -k1,1 -k2,2n -k3,3n | uniq | awk '{print $0"\ter_wMotif`"NR"\t.\t."}'  > ${1%.bed}_peak.bed

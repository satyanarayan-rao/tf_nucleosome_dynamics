#!/bin/bash 
# sh $NGS_SCRIPTS_DIR/run_fimo.sh <fasta> <motif> <out_dir> 
# The only reason I am writing this code is document what fimo command I ran 

fasta=$1
motif=$2
directory=$3
rename=$4

fimo --max-stored-scores 10000000 --oc ${directory} ${motif} ${fasta} 

echo -e "# FIMO command" >> ${directory}/README.md
echo -e "\`\`\`\nfimo --max-stored-scores 10000000 --oc ${directory} ${motif} ${fasta}\n\`\`\`\n" >> ${directory}/README.md
mv ${directory}/fimo.txt $4



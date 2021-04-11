suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("IRanges"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("yaml"))
suppressPackageStartupMessages(library("stringr")) 
suppressPackageStartupMessages(library("R.utils")) 
options(show.error.locations = T)
#options (warn = 2)
source("scripts/utils.R")
args = commandArgs (trailingOnly = T)
# args[1] = tfbs_bed
# args[2] = tcga_count_matrix
# args[3] = output tsv file
# args[4] = chunk size 

# will make the interval tree for the count matrix as it is a 500 bp widow size, and will check for where tfbs belong
tfbs_bed = args[1]
count_matrix_file = args[2]
out_file = args[3]
chunk_size = as.integer (args[4])
count_data = read.csv(count_matrix_file, 
                     sep = "\t", header = FALSE, 
                     stringsAsFactors=FALSE, comment.char="#")


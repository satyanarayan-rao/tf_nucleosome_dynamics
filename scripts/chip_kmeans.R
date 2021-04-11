library(data.table)
library(dplyr)
library(R.utils)
library(ggplot2)
library(ggthemes)
library(stringr)
library(reshape) 
# args[1]: combined data scaled 
# args[2]: out_kmeans_file
# args[3]: number of clusters
# args[4]: output row_order file 
# args[5]: output kmeans center file
# args[6]: output file to store Rdata
options(error=traceback)
args = commandArgs(trailingOnly = T)

merged_df = read.table(args[1], header = T, sep = "\t", stringsAsFactors = F)
only_val_df = merged_df
only_val_df$chrom_loc = NULL
cl = kmeans(only_val_df, centers = as.integer(args[3]), 
            iter.max = 1000, nstart = 20)

merged_df["cl_id"] = cl$cluster 
sorted_merged_df = merged_df[order(merged_df$cl_id), ]
sorted_merged_df$chrom_loc = paste(sorted_merged_df$chrom_loc,
                                   sorted_merged_df$cl_id, sep = "^")
#attach cluster number in 
write.table(sorted_merged_df[,!grepl("cl_id", names(sorted_merged_df))], 
            args[2], row.names = F, col.names = T, sep = "\t", quote = F) 
write.table(sorted_merged_df[, c("chrom_loc", "cl_id")], 
            args[4], row.names = F, col.names = F, sep = "\t", quote = F) 
# write centers in a file
center_df = data.frame(cl$centers)
#center_df$cl_id = row.names(center_df)
write.table(center_df, 
            args[5], sep = "\t", row.names = F, col.names = T, quote = F)
save (cl, file = args[6])


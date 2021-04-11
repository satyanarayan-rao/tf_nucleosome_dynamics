library(data.table)
library(dplyr)
library(R.utils)
library(stringr)
library(reshape) 
# args[1]: input file list 
# args[2]: label list for input files
# args[3]: output file
options(error=traceback)
args = commandArgs(trailingOnly = T)
file_list = unlist(strsplit(args[1], split = " "))
label_list = unlist(strsplit(args[2], split = "@"))
cnt = 1 
df_list = list()
for (f in file_list){
    tmp = read.table(f, stringsAsFactors = F, header = F, sep = "")
    names(tmp) = c ("chrom_loc", label_list[cnt])
    df_list[[label_list[cnt]]] = tmp
    cnt = cnt + 1
}
merged_df = reshape::merge_recurse(df_list, by = "chrom_loc")
to_scale = merged_df[, !grepl("chrom_loc", names(merged_df))] 
scale_df = data.frame(scale(to_scale))
scale_df = cbind(merged_df[, "chrom_loc", drop = F], scale_df) 


write.table(merged_df, 
            args[3], sep = "\t", row.names = F, col.names = T, quote = F)
write.table(scale_df,
            args[4], sep = "\t", row.names = F, col.names = T, quote = F)

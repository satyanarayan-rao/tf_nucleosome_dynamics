library(data.table)
library(dplyr)
library(R.utils)
library(ggplot2)
library(ggthemes)
library(stringr)
library(reshape2) 
library(Cairo)
# args[1]: combined data  
# args[2]: hist pdf file
options(error=traceback)
args = commandArgs(trailingOnly = T)
dt = read.table(args[1], sep = "\t", header = T, stringsAsFactors = F)
dt_sub  = dt[, !grepl("chrom_loc", names(dt))]

dt_sub["to_melt"] = seq(dim(dt_sub)[1]) 
to_plot_df = melt(dt_sub, id.vars = "to_melt")
print (head(to_plot_df))
pdf(args[2]) 
plt = ggplot(to_plot_df, aes(x = value, fill = variable)) + 
      geom_histogram(alpha = 0.5, position = "identity", bins = 50)  + 
      geom_rangeframe() + theme_few()
print (plt)
dev.off()
Cairo::CairoPNG(args[3], height = 4, width = 6, units = "in", res = 150)
print(plt)
dev.off()

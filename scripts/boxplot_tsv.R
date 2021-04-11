library(ggplot2)
library(Cairo)
args = commandArgs (trailingOnly = T)
dt = read.table(args[1], sep = "\t", header = F, stringsAsFactors = F)
dt_sub = dt[, 2, drop = F]
names (dt_sub) = c("Mean.ChIP")
dt_sub$Mean.ChIP = log2(dt_sub$Mean.ChIP + 0.01)
get_cl_id = function(s){return (tail (as.character(unlist(strsplit(s, split="^"))), 1))}
cl_ids = apply(dt[, 1, drop = F], 1, get_cl_id) 
to_plot = cbind(cl_id = cl_ids, dt_sub)
Cairo::CairoPNG(args[2], height = 6, width = 6, res = 150, units = "in")
ggplot (to_plot, aes(x = cl_id, y = Mean.ChIP, fill = cl_id)) + geom_boxplot() + theme_classic() 

dev.off()


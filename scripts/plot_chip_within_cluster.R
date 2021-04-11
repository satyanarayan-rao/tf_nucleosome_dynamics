library(ggplot2)
library(ggthemes)
library(Cairo)
library(stringr)
args = commandArgs(trailingOnly = T)
# args[1]: chip val tsv file
# args[2]: output png 
# args[3]: output pdf
# args[4]: y title

dt = read.table(args[1], sep = "", header = F, stringsAsFactors = F)
names(dt) = c("chr_details", "Mean.ChIP")

dt["pseudo_chip"] = log2 (dt$`Mean.ChIP` + 0.01)

get_cl_id = function(x) {return (tail(unlist(strsplit(x, split = "`")),1) )}
dt["cl_id"] = apply (dt[,"chr_details", drop = F], 1, get_cl_id)

Cairo::CairoPNG(args[2], height = 4, width = 6, units = "in", res = 150)
plt = ggplot(dt, aes(x = cl_id, y = pseudo_chip)) + geom_boxplot(notch = T) + 
      geom_rangeframe() + theme_few() + ylab(args[4])
print(plt)
dev.off()
pdf(args[3], height = 4, width = 6)
print (plt)
dev.off()


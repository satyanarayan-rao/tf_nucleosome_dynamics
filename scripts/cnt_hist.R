library(ggplot2)
library(grid)
library(yaml)
library(ggthemes)
library(scales)
library(stringr)
args = commandArgs(trailingOnly = T) 
options(bitmapType='cairo')

dt = read.table(args[1], sep = "", header = F, stringsAsFactor = F)
names(dt) = c("chrom_details", "cnt")
dt$cnt = as.numeric(dt$cnt)
positive = c("ENCODE", "Carroll", "positive")
negative = c("neg_ctl", "negative")
get_class = function (s) {
              s1 =  str_replace(s, "`neg_ctl", "%negative")
              return (tail(unlist(strsplit(s1, split = "%")),1))}
assign_pos_neg = function (x) {if (x %in% positive) {return ("Pos")} else if (x %in% negative) { return ("Neg")} } 
cl_and_source = apply(dt[, "chrom_details", drop = F], 1, get_class)

dt["cls"] = cl_and_source 

print (head(dt))
png (args[3], height = 4, width = 6, units = "in", res = 150)
plt = ggplot(dt, aes (x = cnt, colour = cls)) + 
      stat_ecdf() + 
      scale_colour_manual(values=c("#69b3a2", "#404080")) + geom_rangeframe() + 
      theme_few() + ggtitle (args[4]) + theme(plot.title = element_text(hjust=0.5)) + 
      xlab("Number of fragments in TFBS") + ylab("eCDF") + xlim(c(0,50))
print (plt)
dev.off()


png (args[2], height = 4, width = 6, units = "in", res = 150)

plt = ggplot(dt, aes (x = cnt, fill = cls)) + 
      geom_histogram(color="#e9ecef", alpha=0.5, position = 'identity', bins = 20) + 
      scale_fill_manual(values=c("#69b3a2", "#404080")) + geom_rangeframe() + 
      theme_few() + ggtitle (args[4]) + theme(plot.title = element_text(hjust=0.5)) + 
      xlab("Number of fragments in TFBS") + ylab("Count") + xlim(c(0,50)) 
print (plt)
dev.off()

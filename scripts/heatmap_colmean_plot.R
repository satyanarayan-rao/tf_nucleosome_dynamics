library(ggplot2)
library(ggthemes)
library(Cairo)
library(yaml)
library(zoo)
args = commandArgs(trailingOnly = T)

dt = read.table(args[1], sep = "\t", header = F, stringsAsFactors = F)

nclust = as.integer(args[5])
colors = yaml.load_file(args[4])
cl_colors = c() 
for (i in seq(nclust)){
    color_reassign =  unlist(colors[i])
    cl_colors = c(cl_colors , color_reassign)
}
print (head(dt))
print (names(cl_colors))
names(cl_colors) = paste0("Cl", names(cl_colors))
####################
figure2_theme <- function (){
    theme (axis.text.y =element_text(vjust =1))+
    theme(plot.title=element_text( size=25 )) +
    theme(axis.title.x = element_text(colour = "black", size = 25),
          axis.title.y = element_text(colour = "black", size = 25)) +
    theme(axis.text.x = element_text(size = 25),
          axis.text.y = element_text(size = 25, vjust = 0.5, hjust = 0.5)) +
    theme(legend.title= element_text(size = 25),
          legend.text = element_text(size = 25)) +
    theme(axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5))
}



dt$rmean = rollmean(dt$V2, k=50, fill = NA)
plt = ggplot(dt, aes(x = V1, y = rmean, color = V3)) + geom_line(lwd = 0.75) + 
      geom_rangeframe() + theme_few() + 
      scale_color_manual (values = cl_colors) + 
      theme(legend.position = "bottom") +  
      figure2_theme() + xlab("Distance from peak summit [bp]") +
      ylab("Enrichment over mean") + 
      guides(color = guide_legend(title = element_blank(), ncol = 1, byrow = F)) + 
      theme(legend.position = c(0.75, 0.6), legend.key = element_rect (fill = "transparent", colour = "transparent"),
            legend.background = element_rect(fill = "transparent", colour = "transparent"))


pdf (args[2], width = 6, height = 6)
print(plt)
dev.off()

Cairo::CairoPNG(args[3], width = 6, height = 6, units = "in", res = 150)
print (plt)
dev.off()


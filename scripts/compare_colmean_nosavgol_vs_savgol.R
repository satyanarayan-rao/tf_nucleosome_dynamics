library(ggplot2)
library(ggthemes)
library(Cairo)
library(yaml)
library(zoo)
args = commandArgs(trailingOnly = T)

dt_nosavgol = read.table(args[1], sep = "\t", header = F, stringsAsFactors = F)
dt_savgol = read.table(args[2], sep = "\t", header = F, stringsAsFactors = F) 

dt_nosavgol$type = "woSavgol" 
dt_savgol$type = "wSavgol" 

nclust = as.integer(args[6])
colors = yaml.load_file(args[5])
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



dt_nosavgol$rmean = dt_nosavgol$V2  #rollmean(dt_nosavgol$V2, k=20, fill = NA)
dt_savgol$rmean = dt_nosavgol$V2 # rollmean(dt_savgol$V2, k=50, fill = NA) # dt_nosavgol$V2 # 

total_data = rbind(dt_nosavgol, dt_savgol)

plt = ggplot(total_data, aes(x = V1, y = rmean, color = V3)) + geom_line(lwd = 0.75) + 
      geom_rangeframe() + theme_few() + 
      scale_color_manual (values = cl_colors) + 
      theme(legend.position = "bottom") +  
      figure2_theme() + xlab("Distance from peak summit [bp]") +
      ylab("Enrichment over mean") + ggtitle(args[7]) +  
      guides(color = guide_legend(title = element_blank(), nrow = 1, byrow = F)) + 
      theme(legend.position = "bottom", legend.key = element_rect (fill = "transparent", colour = "transparent"),
            legend.background = element_rect(fill = "transparent", colour = "transparent"),
           plot.title = element_text(hjust = 0.5))
plt = plt + facet_grid(cols = vars(type))

pdf (args[3], width = 12, height = 6)
print(plt)
dev.off()

Cairo::CairoPNG(args[4], width = 12, height = 6, units = "in", res = 150)
print (plt)
dev.off()


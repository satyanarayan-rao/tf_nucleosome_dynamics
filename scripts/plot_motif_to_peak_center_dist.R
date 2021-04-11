library(ggplot2)
library(ggthemes)
library(cowplot)
library(Cairo)
figure2_theme <- function (){
    theme (axis.text.y =element_text(vjust =1))+
    theme(plot.title=element_text( size=20 )) +
    theme(axis.title.x = element_text(colour = "black", size = 20),
          axis.title.y = element_text(colour = "black", size = 20)) +
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15)) +
    theme(legend.title= element_text(size = 20),
          legend.text = element_text(size = 20)) +
    theme(axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5))
}
args = commandArgs(trailingOnly = T)

dt = read.table(args[1], sep = "\t", header = T, stringsAsFactors = F)

#peak_to_motif_distance	slop	strand
#-113	flank_cnr_200	+
#-74	flank_cnr_200	+
#-10	flank_cnr_200	+
#-93	flank_cnr_200	-
#-104	flank_cnr_200	+
#-108	flank_cnr_200	-
#-74	flank_cnr_200	+
#-70	flank_cnr_200	-
#73	flank_cnr_200	+

dt$abs_val = abs(dt$peak_to_motif_distance)
strand_aware_distance = c() 
cnt_plus_strand = 0 
cnt_minus_strand = 0
for (i in seq (nrow(dt))){
    strand = dt [i, "strand"]
    dist = dt[i, "peak_to_motif_distance"] 
    if (identical(strand, "+")){
        strand_aware_distance = c(strand_aware_distance, dist)
        cnt_plus_strand = cnt_plus_strand + 1 
    }else{
        strand_aware_distance = c(strand_aware_distance, -1*dist)
        cnt_minus_strand = cnt_minus_strand + 1 
    }
}
dt$strand_aware_distance = strand_aware_distance
print(head(dt))

plt_abs = ggplot(dt, aes(x = abs_val)) + 
          geom_histogram( position = "identity", alpha = 0.5, binwidth = 5) + 
          geom_rangeframe() + theme_few() + xlab("Distance from peak summit [bp]") + 
          figure2_theme()  
plt_noabs = ggplot(dt, aes(x = peak_to_motif_distance)) + 
          geom_histogram( position = "identity", alpha = 0.5, binwidth = 5) + 
          geom_rangeframe() + theme_few() + xlab("Distance from peak summit [bp]") + 
          labs(fill = "") + figure2_theme() 

plt_plus_and_minus_strand =  ggplot(dt, aes(x = abs_val, fill = strand)) + 
          geom_histogram( position = "identity", alpha = 0.5, binwidth = 5) + 
          geom_rangeframe() + theme_few() + xlab("Distance from peak summit [bp]") + 
          theme(legend.position = c(0.8, 0.8)) + figure2_theme()

plt_strand_aware_distance = ggplot(dt, aes(x = strand_aware_distance, fill = strand)) + 
          geom_histogram( position = "identity", alpha = 0.5, binwidth = 5) + 
          geom_rangeframe() + theme_few() + xlab("Strand aware up/down-stream [bp]") + 
          theme(legend.position = c(0.8, 0.8)) + figure2_theme()

plt_strand_aware_distance_strands_merged = ggplot(dt, aes(x = strand_aware_distance)) + 
          geom_histogram( position = "identity", alpha = 0.5, binwidth = 5) + 
          geom_rangeframe() + theme_few() + xlab("Strand aware up/down-stream [bp]") + 
          theme(legend.position = c(0.8, 0.8)) + figure2_theme()

total = cnt_plus_strand + cnt_minus_strand 
print(total)
title_lab = paste0 (args[4], " ", args[5], " slop = ", args[6], 
                    "\n(n = ", total, "; + = ", cnt_plus_strand, ", - = ", cnt_minus_strand, ")")
plt = cowplot::plot_grid(plt_abs, plt_strand_aware_distance,plt_noabs, plt_plus_and_minus_strand,
                         nrow = 2, ncol = 2)
tt = ggdraw() + draw_label( title_lab, x = 0.5, hjust = 0.5, size = 25)
 
pdf(args[2], width = 12, height = 12)
print(plot_grid(tt,plt, nrow = 2, rel_heights = c(0.1, 1)))
dev.off() 
Cairo::CairoPNG(args[3], height = 12, width = 12, units = "in", res = 300)
print(plot_grid(tt,plt, nrow = 2, rel_heights = c(0.1, 1)))
dev.off()



plt1 = cowplot::plot_grid(plt_abs, plt_strand_aware_distance_strands_merged,
                         nrow = 1, ncol = 2) 
pdf (args[7], width = 10, height = 5)
print(plot_grid(tt, plt1, nrow = 2, rel_heights = c(0.2, 1)))

dev.off() 
Cairo::CairoPNG (args[8], width = 10, height = 5, units = "in", res = 300)
print(plot_grid(tt, plt1, nrow = 2, rel_heights = c(0.2, 1)))
dev.off() 

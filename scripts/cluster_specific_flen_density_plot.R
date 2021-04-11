library(ggplot2)
library(ggthemes)
library(stringr)
library(Cairo)
args = commandArgs(trailingOnly = T)
file_list = args[1]
cluster_id = as.integer(args[2])
file_list_label = args[3]

#print (file_list, stderr())
#print (cluster_id, stderr())
#print (file_list_label, stderr())
file_vec = unlist(strsplit(file_list, split = " "))
file_label_vec = unlist (strsplit(file_list_label, split = " "))

cnt = 0
to_plot_df = NULL 
for (i in seq(length(file_vec))){
    fname = file_vec[i]
    label = file_label_vec[i]
    tmp_df = read.table(fname, header = F, sep = "\t", stringsAsFactor = F) 
    sub_df = tmp_df[which(tmp_df$V2 == cluster_id), ]
    sub_df$label = label 
    
    if (cnt == 0 ){
        to_plot_df = sub_df
    } else{
        to_plot_df = rbind (to_plot_df, sub_df)
    }
    cnt = cnt + 1 
}
names(to_plot_df) = c("x_tics", "cl_id", "y", "label")
manual_x_tics = seq (35,250, 10)
Cairo::CairoPNG(args[4], height = 4, width = 6, res = 150, units = "in")
ggplot (to_plot_df, aes(x = x_tics, y = y, color = label)) + geom_line()  + geom_rangeframe() + theme_few() +   scale_x_continuous("Fragment length [bp]", breaks = manual_x_tics, labels = as.character(manual_x_tics)) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) + ggtitle(paste0(args[5], "-", "Cluster ", cluster_id)) + theme(plot.title = element_text(hjust = 0.5)) + ylab("Density [A.U.]")
dev.off()

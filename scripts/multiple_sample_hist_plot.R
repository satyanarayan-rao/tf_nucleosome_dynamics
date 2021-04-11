library(ggplot2)
library(ggplot2)
library(ggthemes)
library(stringr)
library(Cairo)
args = commandArgs(trailingOnly = T)
hist_files = args[1]
count_files = args[2]
file_list_label = args[3]

hist_file_vec = unlist(strsplit(hist_files, split = " "))
count_file_vec = unlist(strsplit(count_files, split = " "))
file_label_vec = unlist (strsplit(file_list_label, split = " "))

cnt = 0 
to_plot_df = NULL
for (i in seq(length(hist_file_vec))){
    fname = hist_file_vec[i]
    label = file_label_vec[i]
    tmp_df = read.table(fname, header = F, sep = "\t", stringsAsFactor = F) 
    sub_df = tmp_df[which(tmp_df$V1 >=30 & tmp_df$V1<=250), ]
    
    sub_df$label = label  
    
    if (cnt == 0){
        to_plot_df = sub_df
    } else {
        to_plot_df = rbind(to_plot_df, sub_df)
    }
    cnt = cnt + 1 
    
}
names(to_plot_df) =  c("fragment_length", "hist_val", "sample")
manual_x_tics = seq(30, 250, 10)
Cairo::CairoPNG(args[4], height = 4, width = 6, res = 150, units = "in")
ggplot(to_plot_df, aes(x = fragment_length, y = hist_val, color = sample)) + 
geom_line() + geom_rangeframe() + theme_few() +  
scale_x_continuous("Fragment length [bp]", breaks = manual_x_tics, 
                    labels = as.character(manual_x_tics)) + 
 theme(axis.text.x = element_text(angle = 45, hjust = 1), 
       legend.title = element_blank())
dev.off()


library(ggplot2)
library(ggthemes)
library(reshape2)
options(bitmapType='cairo')
library(scales)
library(yaml)
library(stringr)
source("/beevol/home/satyanarr/workplace/softwares_tools_and_libraries/ngs_scripts/figure_themes.R")
args = commandArgs(trailingOnly = T)
# args[1]: starting mean tsv file that we get from doing mean of k-means kde
# args[2]: manual assignment of cluster id
# args[3]: output meanplot
# args[4]: figure title
# args[5]: gnuplot color file
# args[6]: output manual mean tsv file: just for book-keeping I keep this
# args[7]: peak location file
# args[8]: line plot pdf file
# args[9]: expected length file 
font_size = 25

figure2_theme <- function (){
    theme (axis.text.y =element_text(vjust =1))+
    theme(plot.title=element_text( size=font_size )) +
    theme(axis.title.x = element_text(colour = "black", size = font_size),
          axis.title.y = element_text(colour = "black", size = font_size)) +
    theme(axis.text.x = element_text(size = font_size),
          axis.text.y = element_text(size = font_size)) +
    theme(legend.title= element_text(size = font_size),
          legend.text = element_text(size = font_size)) +
    theme(axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5))
}

mean_tsv_data = read.table(args[1], sep = "\t", 
                  header = F, stringsAsFactors = F, comment.char = "#")
names(mean_tsv_data) = c("x_tics", "cl_id", "mean_density")
manual_cl_id_data = read.table(args[2], header = T, sep ="\t", 
                      stringsAsFactors = F, comment.char = "#")
assigned_cl_id = unlist (lapply(mean_tsv_data$cl_id, function(x){ 
                   manual_cl_id_data[which(manual_cl_id_data$actual == x), "assigned"]}))
# get cl ids and their colors
mean_tsv_data$assigned_cl_id = assigned_cl_id
u_cl_id = sort(unique(mean_tsv_data$assigned_cl_id))
color_dict = yaml.load_file(args[5])
color_list = sapply(u_cl_id, function(x) {color_dict[[as.character(x)]]}) 
print (color_list)
chr_u_cl_id = as.character (u_cl_id)
# read expected length file: 
expt_len_df = read.table(args[9], header =T, sep = "\t") 
expt_len_df$rounded_expt_val = round(expt_len_df$expt)
label_chr_u_cl_id = paste0(chr_u_cl_id, " (", expt_len_df$rounded_expt_val, ")") 
print (chr_u_cl_id)

mean_tsv_data$assigned_cl_id = as.character(mean_tsv_data$assigned_cl_id)
height = 5
width = 6
#print (head(mean_tsv_data))
#print (str(mean_tsv_data))
png (args[3], height = height, width = width, units = "in", res = 300)
kk = ggplot (mean_tsv_data, 
      aes (y = mean_density, x = x_tics, color = assigned_cl_id)) + 
     geom_line(size = 2) + ylab("Mean frequency") + xlab("Fragment length [bp]") + 
     #scale_y_continuous(breaks = round(extended_range_breaks(n=5)(mean_tsv_data$mean_density),4)) + 
     #scale_x_continuous(breaks = extended_range_breaks(n=6)(mean_tsv_data$x_tics)) + 
     geom_rangeframe() + theme_few() + 
     ggtitle(args[4]) + theme(plot.title = element_text(hjust = 0.5)) + 
     scale_color_manual("Cluster\n(W.L.[bp])", breaks = chr_u_cl_id, values = color_list, labels = label_chr_u_cl_id) + 
     figure2_theme() + 
     #theme(legend.position = c(0.85,0.6), legend.key = element_rect (fill = "transparent", colour = "transparent"),
     #       legend.background = element_rect(fill = "transparent", colour = "transparent") ) + 
     theme(legend.position = "bottom", legend.key = element_rect (fill = "transparent", colour = "transparent"),
            legend.background = element_rect(fill = "transparent", colour = "transparent") ) + 
     theme(axis.text.y = element_text(hjust = 0.5, vjust = 0.5)) + 
     guides(color = guide_legend(nrow = 2, byrow = T))

print(kk)
dev.off()
pdf(args[8], height = height, width = width)
print(kk)
dev.off()

# save the manually assigned cluster id dataframe
# get max locations in each cluster
u_assigned_cl_id = str_sort (unique(mean_tsv_data$assigned_cl_id), numeric = T)
locations_of_max = lapply (u_assigned_cl_id, function(x){
                    mean_tsv_data[which.max(
                      mean_tsv_data[which(
                       mean_tsv_data$assigned_cl_id == x), 
                       "mean_density"]), "x_tics"]})
names(locations_of_max) = u_assigned_cl_id 

max_location_df = t(as.data.frame(locations_of_max))
row.names (max_location_df) = u_assigned_cl_id

write.table(mean_tsv_data, args[6], sep = "\t", 
            col.names = T, row.names = F, quote = F)

write.table(max_location_df, args[7], sep = "\t",
            quote = F, row.names = T, col.names = F)

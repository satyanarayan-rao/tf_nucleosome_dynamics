library(ggplot2)
library(ggthemes)
library(reshape2)
#library(Cairo)
options(bitmapType='cairo')
library(scales)
library(yaml)
source("/beevol/home/satyanarr/workplace/softwares_tools_and_libraries/ngs_scripts/figure_themes.R")
# args[1]: input kmeans file (used for heatmap)
# args[2]: output mean plot (png)
# args[3]: x tics 
# args[4]: output mean plot (pdf)
# args[5]: figure title <play with this in snakemake params block>
# args[6]: gnuplot color file
# args[7]: output mean tsv file

args = commandArgs (trailingOnly = T)
input_df = read.table(args[1], sep = "\t", header =F, stringsAsFactors = F)
x_tics = read.table(args[3], header = F, sep="", stringsAsFactors = F)
color_file = args[6]
color_dict = yaml.load_file(color_file)
get_cl_id = function (vv) {return (as.character(unlist(strsplit(vv, split = "\\^"))[1]))}
cl_id = apply(input_df[,1, drop =F], 1, get_cl_id)
with_cl_id = cbind(cl_id = cl_id, input_df)
u_cl_id = sort(as.integer(unique(with_cl_id$cl_id)))
print (u_cl_id, stderr())
mean_vec_list = lapply (u_cl_id, function (x) {apply(with_cl_id[which (with_cl_id$cl_id == as.character(x)), seq(3,dim(with_cl_id)[2])], 2, mean)})
names(mean_vec_list) = u_cl_id
mean_df = as.data.frame(mean_vec_list)
names(mean_df) = u_cl_id
color_list = sapply(u_cl_id, function(x) {color_dict[[x]]}) 
#mean_df["x"] = seq(dim(mean_df)[1]) 
mean_df["x"] = x_tics$V1 
melted_df = melt(mean_df, id.vars = "x") 
names(melted_df) = c("Fragment_length", "Cluster", "Mean_frequency")
# get location of max in each cluster
clusters = unique(melted_df$Cluster) 
max_loc_clusters = lapply(clusters, 
                   function(x) {
                     melted_df[which.max(melted_df[which(melted_df$Cluster == x), "Mean_frequency"]), "Fragment_length"]
                   })
max_freq_clusters = lapply(clusters, 
                   function(x) {
                     max(melted_df[which(melted_df$Cluster == x), "Mean_frequency"])
                   })
rounded_max_loc_clusters = lapply(clusters, 
                   function(x) {
                     round(melted_df[which.max(melted_df[which(melted_df$Cluster == x), "Mean_frequency"]), "Fragment_length"])
                   })

#print(max_loc_clusters)
names(max_loc_clusters) = clusters
names(rounded_max_loc_clusters) = clusters
names (max_freq_clusters) = clusters
golden_ratio = 1.61803398875
height = 6
width = height * golden_ratio
width = 9
# create xtics and color composition 
# start breaks are: 
x_breaks = extended_range_breaks(n=6)(melted_df$Fragment_length)
x_break_colors = rep("black", length(x_breaks))
for (i in clusters){
    if (rounded_max_loc_clusters[[i]] %in% x_breaks){
        # get the location 
        loc = which (x_breaks == rounded_max_loc_clusters[[i]])
        x_break_colors[loc] = color_list[as.integer(i)]
    }else{ 
        x_breaks = c(x_breaks, rounded_max_loc_clusters[[i]])
        x_break_colors = c (x_break_colors, color_list[as.integer(i)])
    }
}
# sort breaks and their corresponding colors 
sort_order = order(x_breaks) 
x_breaks = x_breaks[sort_order]
x_break_colors = x_break_colors[sort_order]
x_breaks_selected = c(x_breaks[1])
x_break_colors_selected = c(x_break_colors[1])
for (i in seq (2, length(x_breaks))){
    if (x_breaks[i]>=(x_breaks_selected[length(x_breaks_selected)] + 10)){
        x_breaks_selected = c(x_breaks_selected, x_breaks[i])
        x_break_colors_selected = c(x_break_colors_selected, x_break_colors[i])
    }else{
        next
    }
} 
print (x_breaks_selected)
print (x_break_colors_selected)

#Cairo::CairoPNG(args[2], height = height, width = width , units = "in", res = 150)
png(args[2], height = height, width = width , units = "in", res = 150)

kk = ggplot (melted_df, aes(y = Mean_frequency, x =Fragment_length, color = Cluster)) + geom_line (size = 2) + 
       ylab("Mean frequency") + xlab("Fragment length [bp]") + 
       #geom_rangeframe() + theme_tufte(base_size = 20, base_family = "Arial", ticks = TRUE) + 
       geom_rangeframe() + theme_few() +  ggtitle(args[5]) + theme(plot.title = element_text(hjust = 0.5)) + 
       scale_y_continuous(breaks = round(extended_range_breaks(n=5)(melted_df$Mean_frequency),4)) + 
       #scale_x_continuous(breaks = extended_range_breaks(n=6)(melted_df$Fragment_length)) + figure1_theme() + 
       scale_x_continuous(breaks = x_breaks_selected) +  
       theme(axis.text.x = element_text(color = x_break_colors_selected),
             axis.ticks.x = element_line(color = x_break_colors_selected)) + figure2_theme() + 
       scale_color_manual(breaks = u_cl_id, values = color_list)  
# Add a vertical line to peak
#print (color_list)
for (i in clusters){
    #print (round(max_loc_clusters[[i]]))
    #print (i)
    #kk = kk + geom_vline (xintercept = max_loc_clusters[[i]], linetype = "dotted")
    kk = kk + geom_segment (x = max_loc_clusters[[i]], y = -Inf, 
                            xend = max_loc_clusters[[i]],yend = max_freq_clusters[[i]],
                            linetype = "dotted", colour = "black")
    #kk = kk + annotate ("text", 
    #                    x = round(max_loc_clusters[[i]]), y = 0, 
    #                    label = as.character(round(max_loc_clusters[[i]])), 
    #                    angle = 45, color = color_list[as.integer(i)])
}
print (kk)
dev.off()
pdf(args[4], height = height, width = width)
print (kk)
dev.off()

# Save melted df 
write.table(melted_df, args[7], col.names = F, 
            row.names = F, quote = F, sep = "\t")

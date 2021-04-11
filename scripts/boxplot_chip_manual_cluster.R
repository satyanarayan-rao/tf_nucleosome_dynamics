library(ggplot2)
library(grid)
library(yaml)
library(ggthemes)
library(scales)
args = commandArgs(trailingOnly = T) 
options(bitmapType='cairo')
# Rscript boxplot_based_on_one_column_with_header.R <input_file> <columd id to use for group plot>  <column id to use for boxplot> <output file>  <pseudo chip> <title>  <label prefix> <x label text> <y label text> 
ks_test = function (v1, v2) {
  return (ks.test(v1, v2)$p.val)
}
figure1_theme <- function (){
    theme (axis.text.y =element_text(vjust =1))+
    theme(plot.title=element_text( size=20 )) +
    theme(axis.title.x = element_text(colour = "black", size = 20),
          axis.title.y = element_text(colour = "black", size = 20)) +
    theme(axis.text.x = element_text(colour = "black",size = 20),
          axis.text.y = element_text(colour = "black", size = 20)) +
    theme(legend.title= element_text(size = 20),
          legend.text = element_text(size = 20)) +
    theme(axis.title.x = element_blank())
}
################### arguments ################
dt = read.table(args[1], sep = "", header = T, stringsAsFactor = F)
#print (head(dt))
#print (args)
dt_colnames = names(dt)
print (dim(dt))
col_to_group = as.numeric(args[2])
col_to_plot = as.numeric(args[3])
out_file = args[4]
pseudo_chip = as.numeric(args[5])
label_title = args[6]
label_prefix = args[7]
x_label = args[8]
y_label = args[9]
ref_for_ks_test = args[10]
min_val = as.integer(args[11])
max_val = as.integer (args[12])
col_dotted_line = as.character (args[13])
ref_box_color = as.character (args[14])
color_file = args [15]
plot_pdf_file = args[16]
plot_title = args[17]
manual_cl_id_data = read.table(args[18], sep = "\t", 
                     header = T, stringsAsFactor = F, comment.char = "#") 
###############################################

dt["Mean.ChIP"] = log2(dt$`Mean.ChIP` + pseudo_chip)
print(head(dt))
dt$Cluster = as.character(dt$Cluster)
#print (head(dt,100))
manual_clusters = unlist(lapply(dt$Cluster, function(x) {
                   manual_cl_id_data[which(manual_cl_id_data$actual == as.integer(x)), "assigned"]}))
dt$Cluster = as.character(manual_clusters)
#print(head(dt,100))
print(manual_clusters[1:10])
print (head(dt))
print (dim(dt))
#dt = dt[which(dt$Mean.ChIP >= -2 & dt$Mean.ChIP <= 5), ]
dt = dt[which(dt$Mean.ChIP >= -8), ]
print(unique(dt$Cluster))
#dt = dt[which(dt$Mean.ChIP >= -6), ]

dt_tmp = dt [, c(col_to_group, col_to_plot)]
#dt_sub = dt [which (dt_tmp$Mean.ChIP >= min_val & dt_tmp$Mean.ChIP <=max_val), ]
dt_sub = dt_tmp
dt_sub_ggplot = dt_sub
colnames (dt_sub_ggplot) = c ("col_to_group", "col_to_plot")
#dt_sub_ggplot = dt_sub_ggplot[which(dt_sub_ggplot$col_to_plot> -6), ]  # discarding all the pseudo-chip value
#print (head(dt_sub))

dt_sub_l = lapply(unique(dt_sub[[dt_colnames[col_to_group]]]), 
                  function (x) { 
		      dt_sub[which(dt_sub[[dt_colnames[col_to_group] ]]==x), dt_colnames[col_to_plot] ] })
names(dt_sub_l) = unique(dt_sub[[dt_colnames[col_to_group]]])
###### Calculate p-values using ref_for_ks_test as reference ########### 

all_ids = unique(dt_sub[[dt_colnames[col_to_group]]])
#cat (all_ids)
#cat ("\n")
ref_for_ks_test_v = c (args[10])
median_val_for_dotted_line = median(dt_sub[which(dt_sub$Cluster == ref_for_ks_test), "Mean.ChIP"])
query_set = sort(setdiff(all_ids, ref_for_ks_test_v))
#print(query_set)

ks_p_vals = lapply (query_set, function(x) {
                       ks_test(unname (dt_sub_l[[x]]), unname(dt_sub_l [[ref_for_ks_test]]) )
                    })
names(ks_p_vals) = query_set
#print(query_set)
print (ks_p_vals)
reassign_zeros = lapply (query_set, function (x) { if (ks_p_vals[[x]] == 0 ) {return("p < 2.2e-16") } else {return(paste0("p = ", round (ks_p_vals[[x]], 4))) }})
p_values_besides_n = lapply (query_set, function (x) { if (ks_p_vals[[x]] == 0 ) {return(paste0("p(",x,",", ref_for_ks_test, ") < 2.2e-16")) } else {return(paste0("p(", x, ",", ref_for_ks_test, ") = ", scientific(ks_p_vals[[x]], digits=2))) }})
names(reassign_zeros) = query_set
to_order_list = c(query_set, ref_for_ks_test)
#color_list = read.table(color_file, header = F, stringsAsFactor = F)$V1 
#color_list = c (color_list,  ref_box_color)
color_dict = yaml.load_file(color_file)
color_list = sapply(query_set, function(x) {color_dict[[x]]})
color_list = unname(c (color_list,  ref_box_color))
print (color_list, stderr())

dt_sub_l_ordered = lapply(to_order_list, function (x) {dt_sub_l[[x]]})
names (dt_sub_l_ordered) = to_order_list 
cat ("######################\n")
print (to_order_list)
cat ("%%%%%%%%%%%%%%%%%%%%%%\n")
print (color_list)


#ggplot_df = data.frame(dt_sub_l_ordered)
dt_sub_ggplot$col_to_group = factor(dt_sub_ggplot$col_to_group, 
                               levels = to_order_list, ordered = T)
png(out_file, height = 10.5, width = 7, units = "in", res = 300)
#pdf(paste0(args[1], "-dist.pdf"), height = 6, width = 9)
max_points = unlist(lapply(dt_sub_l_ordered, max))
min_points = unlist(lapply(dt_sub_l_ordered, min))
median_points = unlist(lapply(dt_sub_l_ordered, median))
counts = unlist(lapply(dt_sub_l_ordered, length))
#boxplot (dt_sub_l_ordered, main = label_title, xlab = x_label, ylab = y_label,
#         frame = F, cex = 0.5, notch = T, ylim = c(min(min_points), max(max_points) + 1.5)) 

# get the upper and lower limits of boxplot y-axis, from boxplot stats
bxplt_dt_sub = boxplot(dt_sub_l_ordered, plot = F) 
#write.table(bxplt_dt_sub, paste0(out_file, ".bxplt.tsv"))
#print (bxplt_dt_sub)
min_for_boxplot = min(bxplt_dt_sub$stats[1,]) 
max_for_boxplot = max(bxplt_dt_sub$stats[5,])
chip_ref = dt_sub[which(dt_sub$Cluster == ref_for_ks_test), "Mean.ChIP"] 
chip_ref = chip_ref[which (chip_ref >= min_for_boxplot & chip_ref <= max_for_boxplot) ]
median_val_for_dotted_line = median(dt_sub[which(dt_sub$Cluster == ref_for_ks_test), "Mean.ChIP"])
#median_val_for_dotted_line = median(chip_ref)

kk = ggplot(dt_sub_ggplot, 
     aes (x = col_to_group, y = col_to_plot, fill = col_to_group, alpha = 0.5)) + 
     geom_boxplot(width = 0.6, notch = T, show.legend = F, position=position_dodge(width=0.6)) + 
     xlab (x_label) + ylab (y_label) + 
     #ylim (c (min_for_boxplot, max_for_boxplot)) + 
     geom_rangeframe() + theme_few() + 
     ggtitle(plot_title) + theme(plot.title = element_text(hjust = 0.5)) + 
     #scale_x_discrete (limits = to_order_list) + 
     scale_fill_manual (breaks = to_order_list, values = color_list) + 
     scale_alpha(guide = 'none') + figure1_theme() + 
     scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
    
# n = on besides each box
for ( i in seq(length(dt_sub_l_ordered)) ){
    #grob <- grobTree(textGrob( paste0("n = ", counts[i]), x= i - 0.5,  y=max_points[i]/2, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))
    if (i < length(dt_sub_l_ordered)){ 
        kk = kk + annotate ("text", x = i - 0.5, y =  median_points[i] + 0.25, #max_points[i]/2, 
                      label = paste0("n = ", counts[i], ", ", p_values_besides_n[i]),
                      angle = 90, size = 5, hjust = 0)}
    else{
        kk = kk + annotate ("text", x = i - 0.5, y = median_points[i] + 0.25, # max_points[i]/2, 
                      label = paste0("n = ", counts[i]), angle = 90, size = 5, hjust=0)
    }
    #kk = kk + annotation_custom(grob)
}
# put p-values and lines 
last = length(dt_sub_l_ordered)
last_minus_one = last - 1
#for (i in seq (last_minus_one)) {
#    kk = kk + geom_segment(x = i, y = max(max_points[i], max_points[last]), 
#             xend = i, yend = max(max_points[i], max_points[last]) + (last - i)*0.5) + 
#    geom_segment(x = last, y = max(max_points[i], max_points[last]), 
#             xend = last, yend = max(max_points[i], max_points[last]) + (last - i)*0.5) + 
#    geom_segment(x = i, y = max(max_points[i], max_points[last]) + (last -i)*0.5, 
#             xend = last, yend = max(max_points[i], max_points[last]) + (last - i)*0.5) + 
#    annotate ("text", (last + i)/2, max(max_points[i], max_points[last])  + (last - i)*0.5 + 0.25, 
#              label = reassign_zeros[[as.character(i)]], size = 4)
#}
kk = kk + geom_hline (yintercept=median_val_for_dotted_line, 
          col = col_dotted_line, lty = 2, size = 1)
print (kk)
#legend ("topright", legend = line_labels_vec, col = color_vals, lty = line_vals)
dev.off()
pdf(args[16], height = 10.5, width = 7)
print (kk)
dev.off()


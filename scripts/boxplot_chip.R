library(ggplot2)
library(grid)
library(yaml)
library(ggthemes)
library(scales)
library(stringr)
library(ggpubr)
source("scripts/expected_fragment_len.R")
args = commandArgs(trailingOnly = T) 




options(bitmapType='cairo')
# Rscript boxplot_based_on_one_column_with_header.R <input_file> <columd id to use for group plot>  <column id to use for boxplot> <output file>  <pseudo chip> <title>  <label prefix> <x label text> <y label text> 
ks_test = function (v1, v2) {
  # v2 is cluster 6 al
  return (ks.test(v2, v1, alternative = "greater")$p.val)
}
figure1_theme <- function (){
    theme (axis.text.y =element_text(vjust =1))+
    theme(plot.title=element_text( size=25 )) +
    theme(axis.title.x = element_text(colour = "black", size = 25),
          axis.title.y = element_text(colour = "black", size = 25)) +
    theme(axis.text.x = element_text(colour = "black",size = 25, hjust=0.5,vjust=0.5),
          axis.text.y = element_text(colour = "black", size = 25, hjust=0.5,vjust=0.5)) +
    theme(legend.title= element_text(size = 25),
          legend.text = element_text(size = 25)) +
    theme(axis.title.x = element_blank())
}


assing_significance_label <- function (p_value) {
    p_values = c (p_value)
    significance_label = c()
    for (p_value in p_values){
        if (p_value> 0.05){
            significance_label = c(significance_label, "n.s.")
        }else if (p_value <=0.05 && p_value > 0.01){
            significance_label = c(significance_label, "*")
        }else if (p_value <=0.01 && p_value > 0.001){
            significance_label = c(significance_label, "**")
        }else if (p_value <=0.001 && p_value > 0.0001){
            significance_label = c(significance_label, "***")
        }else {
            significance_label = c(significance_label, "****")
        }
    }
    return (significance_label)

}

################### arguments ################
dt = read.table(args[1], sep = "", header = T, stringsAsFactor = F)
#print (head(dt))
#print (args)
dt_colnames = names(dt)
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
exact_signal_file = args[18]
log_flag = as.integer(args[24])
median_ratio_plot = args[25]
plt_width = as.integer(args[26])
plt_height = as.integer(args[27])
###############################################
exact_signal_dt = read.table(exact_signal_file, sep ="\t", header = F, stringsAsFactor = F)
names(exact_signal_dt) = c("chrom_details", "Mean.ChIP")  
if (log_flag == 1){ 
    dt["Mean.ChIP"] = log2(dt$`Mean.ChIP` + pseudo_chip) 
    exact_signal_dt["Mean.ChIP"] = log2 (exact_signal_dt$`Mean.ChIP` + pseudo_chip)
}
dt$Cluster = as.character(dt$Cluster)
print (dim (dt))
dt = dt[which(dt$Mean.ChIP > -8), ] 
exact_signal_dt = exact_signal_dt[which(exact_signal_dt$Mean.ChIP > -8), ]
positive = c("ENCODE", "Carroll", "positive")
negative = c("neg_ctl", "negative")
#print(head(exact_signal_dt))
dt["chrom_details"] = exact_signal_dt$chrom_details 
print (head (dt))


get_class = function (s) {
              s1 =  str_replace(s, "`neg_ctl", "%negative")
              return (tail(unlist(strsplit(s1, split = "%")),1))}
assign_pos_neg = function (x) {if (x %in% positive) {return ("Pos")} else if (x %in% negative) { return ("Neg")} } 
cl_and_source = apply(exact_signal_dt[, "chrom_details", drop = F], 1, get_class)
cl_and_source = strsplit(cl_and_source, split = "\\^")
src = unlist(lapply(cl_and_source, function(x) {return (x[1])}))
cl = unlist(lapply(cl_and_source, function(x) {return (x[2])}))
cl_and_source = data.frame(cl = cl, src = src)


pos_neg = lapply(cl_and_source$src, function(x) {assign_pos_neg(x)})
cl_and_source$pos_neg =  pos_neg
cat ("-----------------------------\n")
print(head(cl_and_source))
cat ("-----------------------------\n")
cl_and_source$cl_pos_neg = do.call(paste, c(cl_and_source[, c("cl", "pos_neg")], sep = "-")) 
table_cl_and_source = table(cl_and_source$cl_pos_neg)
ratio_vec = c()
for (i in unique(cl_and_source$cl)){
    pos_count = unname(table_cl_and_source[paste0(i, "-", "Pos")])
    neg_count = unname(table_cl_and_source[paste0(i, "-", "Neg")])
    ratio_vec = c(ratio_vec, round (pos_count/neg_count, 2)) 
}
ratio_df  = data.frame (ratio = ratio_vec, cl = unique(cl_and_source$cl))
print(ratio_df)

#print(table_cl_and_source)

dt_tmp = dt [, c(col_to_group, col_to_plot, 3)]
#dt_sub = dt [which (dt_tmp$Mean.ChIP >= min_val & dt_tmp$Mean.ChIP <=max_val), ]
dt_sub = dt_tmp
dt_sub_ggplot = dt_sub
colnames (dt_sub_ggplot) = c ("col_to_group", "col_to_plot", "chrom_details")
dt_colnames = c ("col_to_group", "col_to_plot", "chrom_details")
#dt_sub_ggplot = dt_sub_ggplot[which(dt_sub_ggplot$col_to_plot> -8), ]  # discarding all the pseudo-chip value
print(head(dt_sub_ggplot))
#dt_sub_ggplot = dt_sub_ggplot[grepl("positive", dt_sub_ggplot$chrom_details), ]
dt_sub = dt_sub_ggplot 
#print (head(dt_sub))

dt_sub_l = lapply(unique(dt_sub[[dt_colnames[col_to_group]]]), 
                  function (x) { 
		                  dt_sub[which(dt_sub[[dt_colnames[col_to_group] ]]==x), 
                      dt_colnames[col_to_plot] ] })
names(dt_sub_l) = unique(dt_sub[[dt_colnames[col_to_group]]])
###### Calculate p-values using ref_for_ks_test as reference ########### 

all_ids = unique(dt_sub[[dt_colnames[col_to_group]]])
print (all_ids)
#cat (all_ids)
#cat ("\n")
ref_for_ks_test_v = c (args[10])
query_set = sort(setdiff(all_ids, ref_for_ks_test_v))
#print(query_set)

ks_p_vals = lapply (query_set, function(x) {
                       ks_test(unname (dt_sub_l[[x]]), 
                       unname(dt_sub_l[[ref_for_ks_test]]) )
                    })
names(ks_p_vals) = query_set
#print(query_set)
print (ks_p_vals)
reassign_zeros = lapply (query_set, function (x) { 
                     if (ks_p_vals[[x]] == 0 ) {return("p < 2.2e-16") } 
                     else {return(paste0("p = ", round (ks_p_vals[[x]], 4))) }})
p_values_besides_n = lapply (query_set, function (x) {
                         if (ks_p_vals[[x]] == 0 ) {
                           return(paste0("p(",x,",", ref_for_ks_test, ") < 2.2e-16")) } 
                         else {return(paste0("p(", x, ",", ref_for_ks_test, ") = ", 
                                         scientific(ks_p_vals[[x]], digits=2))) }})
names(reassign_zeros) = query_set
to_order_list = c(query_set, ref_for_ks_test)
#color_list = read.table(color_file, header = F, stringsAsFactor = F)$V1 
#color_list = c (color_list,  ref_box_color)
color_dict = yaml.load_file(color_file)
color_list = sapply(query_set, function(x) {color_dict[[x]]})
color_list = unname(c (color_list,  ref_box_color))
print (color_list, stderr())
print (to_order_list)

dt_sub_l_ordered = lapply(to_order_list, function (x) {dt_sub_l[[x]]})
names (dt_sub_l_ordered) = to_order_list 
cat ("######################\n")
print (to_order_list)
cat ("%%%%%%%%%%%%%%%%%%%%%%\n")
print (color_list)


#ggplot_df = data.frame(dt_sub_l_ordered)
dt_sub_ggplot$col_to_group = factor(dt_sub_ggplot$col_to_group, 
                               levels = to_order_list, ordered = T)
png(out_file, height = plt_height, width = plt_width, units = "in", res = 300)
#pdf(paste0(args[1], "-dist.pdf"), height = 6, width = 9)
max_points = unlist(lapply(dt_sub_l_ordered, max))
min_points = unlist(lapply(dt_sub_l_ordered, min))
median_points = unlist(lapply(dt_sub_l_ordered, median))
print ("############ %%%%%%%%%%%%% ############")
print (median_points)
print(min_points)

######################### use expected fragment length in cluster to assing cluster id #################### 
ret_list = expected_flen(args[21])

order_by_expected_flen = ret_list$e_df
verbose_dt = ret_list$new_dt
write.table(order_by_expected_flen, file = args[22],
            col.names = T, row.names = F, quote = F, sep = "\t")
write.table(verbose_dt, file = args[23],
            col.names = T, row.names = F, quote = F, sep = "\t")



################################# 

#median_df = data.frame(actual = order(median_points, decreasing = T),
#                       assigned = seq(length(median_points)), 
#                       stringsAsFactors = F) 
median_df = data.frame(actual = order_by_expected_flen$actual, 
                       assigned = order_by_expected_flen$assigned,
                       stringsAsFactors = F)
median_df = median_df[order(median_df$actual), ]
ratio_df["assigned"] = unlist (lapply(ratio_df$cl, function(x){
                         return (median_df[median_df$actual == x, "assigned"] )})) 
ratio_df = ratio_df[order(ratio_df$assigned), ]
ratio_df["assigned"] = as.character(ratio_df$assigned)

assigned_cl_id = unlist(lapply (dt_sub_ggplot$col_to_group, function(x){
                    return (median_df[median_df$actual == x, "assigned"])}))
dt_sub_ggplot["assigned_cl_id"] = as.character(assigned_cl_id)
median_val_for_dotted_line = median(dt_sub_ggplot[which(dt_sub_ggplot$assigned_cl_id == ref_for_ks_test), "col_to_plot"])
print(head(dt_sub_ggplot))
dt_sub_ggplot$assigned_cl_id = factor(dt_sub_ggplot$assigned_cl_id, 
                               levels = to_order_list, ordered = T)
ratio_df$assigned = factor(ratio_df$assigned, 
                          levels = to_order_list, ordered = T)
print (str(dt_sub_ggplot))

#dt_sub_ggplot = dt_sub_ggplot[order(dt_sub_ggplot$assigned_cl_id), ]
write.table(median_df, args[20], col.names = T, row.names = F, sep = "\t", quote = F)  
#counts = unlist(lapply(dt_sub_l_ordered, length))
counts_in_group = table(dt_sub_ggplot$assigned_cl_id) 
counts = rep(0, length(counts_in_group))
for (i in names(counts_in_group)){
    counts[as.integer(i)] = unname(counts_in_group[i])
}
# since I have reassigned cluster number - I have to calculated p-value here with reassigned clusters 
#print ("------------SSSSSSSSSSS-----------")
#print(c(names()))
ks_p_vals = lapply (query_set, function(x) {
                       ks_test(dt_sub_ggplot[which(dt_sub_ggplot$assigned_cl_id == as.character(x)), "col_to_plot"], 
                               dt_sub_ggplot[which(dt_sub_ggplot$assigned_cl_id == as.character(ref_for_ks_test)), "col_to_plot" ])
                    }) 
names(ks_p_vals) = query_set 
reassign_zeros = lapply (query_set, function (x) { 
                     if (ks_p_vals[[x]] == 0 ) {return("p < 2.2e-16") } 
                     else {return(paste0("p = ", round (ks_p_vals[[x]], 4))) }})
p_values_besides_n = lapply (query_set, function (x) {
                         if (ks_p_vals[[x]] == 0 ) {
                           return(paste0("p(",x,",", ref_for_ks_test, ") < 2.2e-16")) } 
                         else if (ks_p_vals[[x]] > 0.05) {
                           return ("")}
                         else {return(paste0("p(", x, ",", ref_for_ks_test, ") = ", 
                                         scientific(ks_p_vals[[x]], digits=2))) }}) 
#print ("----------SATYA----------\n")
print (p_values_besides_n)

#boxplot (dt_sub_l_ordered, main = label_title, xlab = x_label, ylab = y_label,
#         frame = F, cex = 0.5, notch = T, ylim = c(min(min_points), max(max_points) + 1.5)) 
#
# get the upper and lower limits of boxplot y-axis, from boxplot stats
bxplt_dt_sub = boxplot(dt_sub_l_ordered, plot = F) 
#write.table(bxplt_dt_sub, paste0(out_file, ".bxplt.tsv"))
#print (bxplt_dt_sub)
min_for_boxplot = min(bxplt_dt_sub$stats[1,]) 
max_for_boxplot = max(bxplt_dt_sub$stats[5,])
#chip_ref = dt_sub[which(dt_sub$Cluster == ref_for_ks_test), "Mean.ChIP"] 
#chip_ref = chip_ref[which (chip_ref >= min_for_boxplot & chip_ref <= max_for_boxplot) ]
#median_val_for_dotted_line = median(dt_sub[which(dt_sub$col_to_group == ref_for_ks_test), "col_to_plot"])
#median_val_for_dotted_line = median(chip_ref)
cat("---------------------------------##\n")
print (head(dt_sub_ggplot)) 
print(names(dt_sub_ggplot))
cat("---------------------------------##\n")

cat ("------------------ medians -------------\n")
print (c(min_for_boxplot, max_for_boxplot))
print (max(dt_sub_ggplot$col_to_plot))
for_median = dt_sub_ggplot[which(dt_sub_ggplot$assigned_cl_id == ref_for_ks_test), "col_to_plot"]
median_val_for_dotted_line = median(for_median[for_median<=max_for_boxplot & for_median >=min_for_boxplot])
kk = ggplot(dt_sub_ggplot, 
     aes (x = assigned_cl_id, y = col_to_plot, fill = assigned_cl_id, alpha = 0.5)) + 
     geom_boxplot(width = 0.6, notch = T, show.legend = F, position=position_dodge(width=0.6)) + 
     xlab (x_label) + ylab (y_label) + 
     #ylim (c (min_for_boxplot, max_for_boxplot)) + 
     geom_rangeframe() + theme_few() + 
     ggtitle(plot_title) + theme(plot.title = element_text(hjust = 0.5)) + 
     #scale_x_discrete (limits = to_order_list) + 
     scale_fill_manual (breaks = to_order_list, values = color_list) + 
     scale_alpha(guide = 'none') + figure1_theme() + 
     scale_y_continuous(labels = scales::number_format(accuracy = 0.1))  # + ylim (c (-4, 4))
     #ylim (c (min_for_boxplot, max_for_boxplot)) 
     #ylim (c (100, 700)) 

if (identical(min_val, as.integer("-1000"))){
    kk = kk + ylim(c(min_for_boxplot, max_for_boxplot)) 
}else{
    print ("%%%%%%%%%%%%%%%%%%--------- Check ---------------- %%%%%%%%%%%")
    kk = kk + ylim(c(min_val, max_val)) 
} 
    
# n = on besides each box
bonferroni_list = c()
bonferroni_sig_str = c()
n_p_text = c()
for ( i in seq(length(dt_sub_l_ordered)) ){
    #grob <- grobTree(textGrob( paste0("n = ", counts[i]), x= i - 0.5,  y=max_points[i]/2, hjust=0, gp=gpar(col="red", fontsize=13, fontface="italic")))
    if (i < length(dt_sub_l_ordered)){ 
        #kk = kk + annotate ("text", x = i - 0.5, y =  min_for_boxplot + 0.25, #max_points[i]/2, 
        #              label = paste0("n = ", counts[i], " ", p_values_besides_n[i]),
        #              angle = 90, size = 5, hjust = 0)  
        if (ks_p_vals[[i]] == 0){
            n_p_text = c(n_p_text, paste0("n = ", counts[i], ", ", "p (", i, ",6) < ", "2.2e-16")) 
        }else{
            n_p_text = c(n_p_text, paste0("n = ", counts[i], ", ", "p (", i, ",6) = ", format(ks_p_vals[[i]], digits = 2) )) 
        }
        # get bonferroni corrected p-values
        if (ks_p_vals[[i]] == 0){
            p_val = 2e-16
            bonferroni_list = c(bonferroni_list, 2e-16)
            s_str = assing_significance_label (p_val)
            kk = kk + annotate ("text", x = i, y =  max_for_boxplot , #max_points[i]/2, 
                      label = s_str, size = 7, hjust = 0.5)  
 
        }else{
            p_val = counts[i]*ks_p_vals[[i]]
            bonferroni_list = c(bonferroni_list, p_val)
            s_str = assing_significance_label (p_val)
            if (s_str != "n.s."){
            kk = kk + annotate ("text", x = i, y =  max_for_boxplot , #max_points[i]/2, 
                      label = s_str, size = 7, hjust = 0.5) 
            }
            
             
        }
        
        pos_count = unname(table_cl_and_source[paste0(i, "-", "Pos")])
        neg_count = unname(table_cl_and_source[paste0(i, "-", "Neg")])
        ratio = round(pos_count/neg_count, 2)
      
   }else{
        #kk = kk + annotate ("text", x = i - 0.5, y = min_for_boxplot + 0.25, # max_points[i]/2, 
        #              label = paste0("n = ", counts[i]), angle = 90, size = 5, hjust=0)
        pos_count = unname(table_cl_and_source[paste0(i, "-", "Pos")])
        neg_count = unname(table_cl_and_source[paste0(i, "-", "Neg")])
        n_p_text = c(n_p_text, paste0("n = ", counts[i]))
    }
    #kk = kk + annotation_custom(grob)
}
bonferroni_sig_str = assing_significance_label(bonferroni_list)
bonferroni_df = data.frame(b_p_val = bonferroni_list, b_s_str = bonferroni_sig_str)
n_p_df = data.frame(np = n_p_text) 
# put p-values and lines 
last = length(dt_sub_l_ordered)
last_minus_one = last - 1
dt_med_ratio_bxplt = NULL
if (identical(min_val, as.integer("-1000"))){

    kk = kk + geom_hline (yintercept=median_val_for_dotted_line, 
              col = col_dotted_line, lty = 2, size = 1) 
    med_ratio_vec = c() 
    cl_id_vec = c()
    for (cl_id in as.character(seq(as.integer(ref_for_ks_test)))){
        med_val = median(dt_sub_ggplot[(dt_sub_ggplot$assigned_cl_id == cl_id) & (dt_sub_ggplot$col_to_plot >= min_for_boxplot) & (dt_sub_ggplot$col_to_plot <= max_for_boxplot), "col_to_plot"]) 
        med_ratio_vec = c(med_ratio_vec, med_val)
        cl_id_vec =  c(cl_id_vec, as.character(cl_id))
    }
    dt_med_ratio_bxplt = data.frame(cl_id = cl_id_vec, med_ratio = med_ratio_vec/med_ratio_vec[as.integer(ref_for_ks_test)], med_val = med_ratio_vec)
}else{
    tmp_med  = median(dt_sub_ggplot[(dt_sub_ggplot$col_to_plot >= min_val) &(dt_sub_ggplot$col_to_plot <= max_val) & (dt_sub_ggplot$assigned_cl_id == ref_for_ks_test), "col_to_plot"])
     kk = kk + geom_hline (yintercept=tmp_med, 
              col = col_dotted_line, lty = 2, size = 1)    
    med_ratio_vec = c()
    cl_id_vec = c()
    for (cl_id in as.character(seq(as.integer(ref_for_ks_test)))){
        med_val = median(dt_sub_ggplot[(dt_sub_ggplot$col_to_plot >= min_val) &(dt_sub_ggplot$col_to_plot <= max_val) & (dt_sub_ggplot$assigned_cl_id == cl_id), "col_to_plot"]) 
        med_ratio_vec = c(med_ratio_vec, med_val)
        cl_id_vec =  c(cl_id_vec, as.character(cl_id))
    }
    dt_med_ratio_bxplt = data.frame(cl_id = cl_id_vec, med_ratio = med_ratio_vec/med_ratio_vec[as.integer(ref_for_ks_test)], med_val = med_ratio_vec)
    
}
print (kk)
dev.off() 
pdf(args[16], height = plt_height, width = plt_width)
print (kk)
dev.off()

write.table (bonferroni_df, paste0(args[16], ".sig.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(n_p_df, paste0(args[16], ".np.tsv"), row.names = F, col.names = F, quote = F) 

print (head (ratio_df))
print (args)
print (args[19])

png(args[19], width = 5, height = 9.5, units = "in", res = 150)
#ratio_df  = data.frame (ratio = ratio_vec, cl = unique(cl_and_source$cl) 
plt = ggplot(ratio_df, aes(x = assigned, y = ratio)) +
       geom_bar(stat = "identity", width = 0.5) + 
       geom_rangeframe() + theme_few()  + figure1_theme() + 
       geom_hline(yintercept = 1, col = "red", lty = 2) + 
       ylab("P/N")
qq = ggarrange(kk, plt, ncol = 1, nrow = 2, heights = c(5,1))
print(qq)
dev.off()


## do the median ratio plot now
print ("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
print (dt_med_ratio_bxplt)
median_ratio_plt = ggplot(dt_med_ratio_bxplt, aes(x = cl_id, y = med_ratio )) +
    geom_bar(stat = "identity") + geom_rangeframe() + theme_few() + figure1_theme()
pdf(args[25], width =6, height = 5)
print(median_ratio_plt)
dev.off()

write.table(dt_med_ratio_bxplt, file = args[28], sep = "\t", 
            row.names = F, col.names = F, quote = F)



print("@@@@@@@@@@@@@@@@@@@@")
print (c(min_for_boxplot, max_for_boxplot, median_val_for_dotted_line)) 

write.table(dt_sub_ggplot, file = "tmp/jnk.tsv", sep = "\t", row.names = F, quote = F, col.names =  F)

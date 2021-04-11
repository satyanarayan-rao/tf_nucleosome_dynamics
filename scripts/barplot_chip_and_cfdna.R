library(data.table)
library(dplyr)
library(R.utils)
library(stringr)
library(reshape) 
library(ggplot2)
library(ggthemes)
library(Cairo)
library(cowplot)
options(error=traceback)

args = commandArgs(trailingOnly = T)
file_list = unlist(strsplit(args[1], split = " "))
label_list = unlist(strsplit(args[2], split = "@"))
print (file_list) 
print(label_list)
cnt = 1 
df_list = list()
df_list_no_name_changed = list()
for (f in file_list){
    tmp = read.table(f, stringsAsFactors = F, header = T, sep = "")
    tmp = tmp [!grepl("not.found", tmp$cfdna_cl_id), ]
    tmp1 = tmp 
    print(table(tmp$chip_kmean_cl_id))
    to_modify_header = names(tmp) 
    for (j in seq(2, length(to_modify_header))){
        to_modify_header[j] = paste(label_list[cnt], 
                               to_modify_header[j], sep = "_")
    }
      
    names(tmp) =  to_modify_header
    df_list[[label_list[cnt]]] = tmp 
    df_list_no_name_changed[[label_list[cnt]]] = tmp1
    cnt = cnt + 1
}
#merged_df = reshape::merge_recurse(df_list, by = "chr_loc")
merged_df =  inner_join(df_list[[label_list[1]]], df_list[[label_list[2]]], by = "chr_loc")
total_inner = dim(merged_df)[1]
write.table(merged_df, "tt.tsv", row.names = F, col.names = T, sep = "\t", quote = F)
print(head(merged_df))


# percentage of sites falling in chip clusters

chip_percentage_cl = (table(merged_df[,paste0(label_list[1], "_chip_kmean_cl_id")])/length(merged_df[,paste0(label_list[1], "_chip_kmean_cl_id")]))*100
chip_count_cl = table(merged_df[,paste0(label_list[1], "_chip_kmean_cl_id")])

table_df_list = list()
for (l in label_list){
    merged_df[paste0(l, "_and_chip")] = paste(merged_df[,paste0(l, "_cfdna_cl_id") ],
                                              merged_df[,paste0(l, "_chip_kmean_cl_id") ],
                                              sep = "-") 
    tmp = data.frame(table(merged_df[paste0(l, "_and_chip")]))
    tmp[l] = unlist(lapply(as.character(tmp$Var1), function(x){
                 unlist(strsplit(x, split = "-"))[1]
             }))
    tmp["chip_cluster"] = unlist(lapply(as.character(tmp$Var1), function(x){
                 unlist(strsplit(x, split = "-"))[2]
             })) 
    
    tmp$rel_frequency = unlist(lapply(seq(dim(tmp)[1]), function (x) {
        rel = tmp[x, "Freq"]/unname(chip_count_cl[tmp[x, "chip_cluster"]])
        return(rel)
       })) 
    table_df_list[[l]] = tmp
}
plt_list = list() 

for (l in label_list){
    plt_list[[l]] = ggplot(table_df_list[[l]], 
                       aes_string(x = l, y = "rel_frequency", 
                                  fill = "chip_cluster" )) + 
          geom_bar(stat = "identity", position = "dodge") + geom_rangeframe() +
          theme_few() + ylab("Fraction of sites") + ggtitle(paste0("Total sites = ", total_inner)) + theme(plot.title = element_text(hjust = 0.5))
    
}
Cairo::CairoPNG(args[3], height = 4, width = 12, units = "in", res = 150)
gl = lapply(plt_list, ggplotGrob)
plot_grid(plotlist = gl, nrow = 1, ncol = 2, align = "h")
dev.off() 


# go solo - as in just take df_list and plot 
nf_removed_list = list()
nf_plt_list = list()
for (l in label_list){ 
    tmp_df = df_list_no_name_changed[[l]]
  
    tmp_df = tmp_df[!grepl("not.found", tmp_df$cfdna_cl_id), ]
    t_sites = dim(tmp_df)[1]
    tmp_df$c_d = paste(tmp_df$chip_kmean_cl_id, tmp_df$cfdna_cl_id, sep = "-")
    print(head (tmp_df))
    chip_count = table(tmp_df$chip_kmean_cl_id)
    kt = data.frame(table(tmp_df$c_d)) 
    kt$chip = unlist(lapply(as.character(kt$Var1), function(x){
      unlist (strsplit(x, split = "-"))[1]
    }))
    kt[l] = unlist(lapply(as.character(kt$Var1), function(x){
      unlist (strsplit(x, split = "-"))[2]
    }))
    kt$rel_f = unlist(lapply( seq(dim(kt)[1]), function(x){
      rel = kt[x, "Freq"]/ chip_count[kt[x, "chip"]]
      return (rel)
    }))
    nf_removed_list[[l]] = kt
    #print(head(kt))
    plt = ggplot(kt, aes_string (x = l, y = "rel_f", fill = "chip")) + 
          geom_bar(stat = "identity", position = "dodge") +  geom_rangeframe() + 
          theme_few() + ylab ("Fraction of sites") + labs(fill = "chip_cluster")  + 
          ggtitle(paste0("Total sites = ", t_sites)) + theme(plot.title = element_text(hjust = 0.5))
    nf_plt_list[[l]] = plt
}

Cairo::CairoPNG(paste0(args[3], "-individual.png"), height = 4, width = 12, units = "in", res = 150)
gl = lapply(nf_plt_list, ggplotGrob)
plot_grid(plotlist = gl, nrow = 1, ncol = 2, align = "h")
dev.off() 

library(data.table)
library(dplyr)
library(R.utils)
library(ggplot2)
library(ggthemes)
options(error=traceback)
args = commandArgs(trailingOnly = T)

dt = read.table(args[1], sep = "", header = F, 
                stringsAsFactors = F)  
row.names(dt) = dt$V1
dt_only_values = dt[, 2:length(dt)]
set.seed (56)
num_cluster = as.integer(args[4])
cl = kmeans(dt_only_values, centers = num_cluster, 
            iter.max = 1000, nstart = 20) 
save(cl, file = args[5])
res_data = dt

#order clusters based on mean around the center of `center vectors` 
# order_of_cluster centers
num_cols_cl_center = dim(cl$centers)[2] 
mean_center_of_centers = NULL 
res_data ["cl_id"] = cl$cluster
ordered_res_data = NULL
#if (num_cols_cl_center/2 - 50 > 0) {
#    left = as.integer(num_cols_cl_center/2) - 50
#    right = as.integer(num_cols_cl_center/2) + 50
#    mean_center_of_centers = order(rowMeans(cl$centers[, c(left:right)]))
#    # we are defining new cluster names based on mean of the centers
#    # Clsuster 1 should be the one with center mean as most negative
#    # characterizing of nucleosome absence.
#    tmp_df = data.frame (mean_center_of_centers)
#    to_order_df = tmp_df[order(tmp_df$mean_center_of_centers), , drop = F]
#    to_order_df["actual_cl_id"] = row.names (to_order_df)
#    to_order_df["assigned_cl_id"] = c (1:dim(to_order_df)[1])
#    assigned_cl_id_vec = unname(sapply(cl$cluster, function(x){ to_order_df[x, "assigned_cl_id"]}))
#    res_data["assigned_cl_id"] = assigned_cl_id_vec
#    res_data["rmean_center_100"] = rowMeans (res_data[, c(left:right)])
#    ordered_res_data = res_data[order(res_data$assigned_cl_id, res_data$rmean_center_100), ]
#}
######### using just the centers ###################################
center_of_centers = data.frame(cl$centers[, as.integer(dim(cl$centers)[2]/2), drop = F ]) 
colnames(center_of_centers) = c("middle")
center_of_centers["actual"] = row.names(center_of_centers) 
reordered_center_of_centers = center_of_centers[order(center_of_centers$middle), , drop = F]
reordered_center_of_centers["value_ordered"] = seq(dim(reordered_center_of_centers)[1]) 
write.table (reordered_center_of_centers, paste0(args[3], "_reoder_center.tsv"), 
             col.names =T, row.names = T, quote = F)

write.table(cl$centers, paste0(args[3], "_cluster_center_vec.tsv"), col.names = F, row.names = F, quote = F)

#center_of_centers["assigned_cl_id"] = as.integer (row.names(reordered_center_of_centers))

assigned_cl_id_vec = unname(sapply(cl$cluster, function(x){ reordered_center_of_centers[which (reordered_center_of_centers$actual == x), "value_ordered"]}))
#print (assigned_cl_id_vec)
res_data["assigned_cl_id"] = assigned_cl_id_vec
left = as.integer(num_cols_cl_center/2) - 25
right = as.integer(num_cols_cl_center/2) + 25
res_data["rmean_center_50"] = rowMeans (res_data[, c(left:right)])
#print (center_of_centers)
#print (res_data[1:2, c(10, 11, 402, 403)])
#print (colnames(res_data))
ordered_res_data = res_data[order(-res_data$assigned_cl_id), ]
#write.table(ordered_res_data, paste0(args[2], "_test.tsv"),
#            col.names = F, row.names = F, sep = "\t", quote = F)
# delete `cl_id` column
ordered_res_data["cl_id"] = NULL

####### Write the rownames in tsv ############# 
tmp_df = data.frame (row_names = row.names (ordered_res_data), cluster_id = ordered_res_data$assigned_cl_id)
write.table (tmp_df, args[3], 
             quote = F, row.names = F, col.names = F, sep = "\t") 
#####################################################################
ordered_res_data["rmean_center_50"] = NULL
ordered_res_data["assigned_cl_id"] = NULL
##### Write gzip file after k-means ########## 
out_to_write = tools::file_path_sans_ext(args[2])
write.table(ordered_res_data, file = out_to_write, 
            col.names = F, row.names = F, sep = "\t", quote = F)
gzip(out_to_write, overwrite = T)  


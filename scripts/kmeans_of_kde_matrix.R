args = commandArgs(trailingOnly = T)
options (traceback = T)

num_cluster = as.numeric (args[2])
max_iter = as.numeric (args[3])
dt = read.table (args[1], sep ="\t", header = F, 
                               stringsAsFactors = F)
row_cols_dt = dim(dt)
dt_sub = dt[, seq(5, row_cols_dt[2])]
print (row_cols_dt)
print (head(dt_sub[1:5, 1:5]))
set.seed(731)
kmeans_dt = kmeans(dt_sub, centers = num_cluster, 
                   iter.max = max_iter,  nstart = 20)
##### Reorder clusters by index of max and max (both are required in case of ties max will be used)
max_of_centers = apply(kmeans_dt$centers, 1, max)
idx_of_max_of_centers = apply(kmeans_dt$centers, 1, which.max)
idx_and_max_df = data.frame(idx = idx_of_max_of_centers, max_vals = max_of_centers)
idx_and_max_df["actual"] = row.names(idx_and_max_df)
ordered_by_idx_and_max = idx_and_max_df[order(idx_and_max_df$idx, idx_and_max_df$max_vals), ]
ordered_by_idx_and_max["assigned_cluster_id"] = seq (dim(ordered_by_idx_and_max)[1]) 
get_assigned_cl_id = function (actual_cl_id, center_df){
    return(center_df[which(center_df$actual == actual_cl_id), "assigned_cluster_id"])
}
assigned_cl_id_vec = unname(sapply(kmeans_dt$cluster, function(x){
                     return(get_assigned_cl_id(x, ordered_by_idx_and_max))}))


out_dt = dt
#out_dt = cbind(cl_id = kmeans_dt$cluster, out_dt)
out_dt = cbind(cl_id = unname(assigned_cl_id_vec), out_dt)

ordered_out_dt = out_dt[order(-out_dt$cl_id), ]
#print (head(ordered_out_dt[1:5,1:5]))
merged_columns = do.call(paste, c(ordered_out_dt[, seq(5)], sep = "^"))
print (merged_columns[1:5])
bed_data_with_cl_id = do.call(paste, c(ordered_out_dt[, seq(5)], sep = "\t"))
# nullify the first five columns, as they have been merged
ordered_out_dt[, seq(5)] = NULL
ordered_out_dt = round (ordered_out_dt, 4)
ordered_out_dt = cbind(cl_id_and_chr = merged_columns, ordered_out_dt)

write.table(ordered_out_dt, file = args[4], row.names = F, 
            col.names = F, sep = "\t", quote = F)
write.table(bed_data_with_cl_id, file = args[5], 
            row.names = F, col.names = F, quote = F, sep = "\t")


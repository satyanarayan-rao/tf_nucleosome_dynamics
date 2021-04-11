library(MASS)
args = commandArgs(trailingOnly = T)
options(traceback =T)
############# Parse arguments ################### 
out_file = args[2]
kde_n = as.numeric(args[3])
kde_bw = as.numeric(args[4])
min_flen = as.numeric(args[5])
max_flen = as.numeric(args[6])
################# read data #####################
dt = read.table (args[1], sep = "", header = F, stringsAsFactors = F)
row_cols_dt = dim (dt)
dt_sub = dt[, seq(5, row_cols_dt[2])]
#### kde2d and flatten it by row (meaning that frequency at flen range) ### 
kde2d_and_flatten = function (flen_count_vec, x_points, x_total_bins, y_total_bins){
    
    kde_data = MASS::kde2d(x = x_points, y = flen_count_vec, h = c(5, 5), n = c(x_total_bins, y_total_bins))
    averge_over_y = round (apply(kde_data$z, 1, sum), 4)
    ret_list = list()
    ret_list["y"] = averge_over_y 
    ret_list["x"] = kde_data$x
    return (ret_list)
} 

#print(dt_sub[1:5, 1:5])
row_wise_kde = apply(dt_sub, 1, kde2d_and_flatten, 
                     x_points = seq(min_flen, max_flen), 
                     x_total_bins = as.integer((max_flen-min_flen)/5),
                     y_total_bins = 5) 
prnt 
names(row_wise_kde) = row.names (dt_sub)
just_one = row.names(dt_sub)[1]
only_density_values = lapply (names(row_wise_kde), function (x) {row_wise_kde[[x]]$y})
names(only_density_values) = names (row_wise_kde)
density_df = t(round (as.data.frame(only_density_values), 4))
names (density_df) = as.character (row_wise_kde[just_one]$x)

# put the bed details in the beginning
density_df = cbind(dt[, seq(1,4)], density_df)
write.table(density_df, file = args[2], row.names = F, col.names = F, quote =F, sep = "\t")

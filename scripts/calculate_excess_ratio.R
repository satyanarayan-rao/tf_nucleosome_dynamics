args = commandArgs(trailingOnly = T)
dt = read.table(args[1], sep = "\t", header = F, stringsAsFactors = F)
cluster_wise_count_a = apply(dt, 1, sum)
cluster_wise_count_b = apply(dt, 2, sum)
total_a = sum(cluster_wise_count_a) 
total_b = sum(cluster_wise_count_b) 
print (c(total_a, total_b))

if (total_a == total_b){
    obs = dt/total_b
    print(obs)
    expected = (cluster_wise_count_a%o%cluster_wise_count_b)/(total_a*total_b) 
    print(expected)
    obs_by_exp = round(obs/expected,2)
    ## write the exess ratio matirx 
    write.table(obs_by_exp, file = args[2], sep = "\t", row.names = F, col.names = F, quote = F)
    
    
}else{
    stop("Sites count should be the same ... first find common TFBSs in both systems! Exiting...")
}

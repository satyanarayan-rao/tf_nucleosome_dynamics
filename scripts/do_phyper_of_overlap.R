library(stats)
args = commandArgs(trailingOnly = T) 
assing_significance_label <- function ( p_values ) {

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
count_dt = as.matrix(read.table(args[1], header = F, stringsAsFactors = F, sep = "\t"))
excess_ratio_dt = as.matrix(read.table(args[4], header = F, stringsAsFactors = F, sep = "\t" ))
rsum = apply(count_dt, 1, sum)
csum = apply(count_dt, 2, sum)

total_c = sum(csum)
total_r = sum(rsum)
print(count_dt)
print (nrow(count_dt))
print (ncol(count_dt))

p_dt = count_dt
s_dt = count_dt
p_vec = c()
for (i in seq(nrow(count_dt))){
    for (j in seq(ncol(count_dt))){
        q = count_dt[i,j]
        m = unname(rsum[i])
        n = total_r 
        k = unname(csum[j])
        pval = NULL
        if (excess_ratio_dt[i, j] >= 1){
            pval = phyper(q, m, n, k, lower.tail = F, log.p = F)
        }else{
            pval = phyper(q, m, n, k, lower.tail = T, log.p = F)
        }
        p_dt[i,j] = pval
        s_dt[i,j] = assing_significance_label(pval)
        p_vec = c(p_vec, pval)
        

    }
}

adj_p = p.adjust(p_vec, method = "fdr") 
adj_p_df = count_dt
adj_s_df = count_dt
for (i in seq(nrow(count_dt))){
    for (j in seq(nrow(count_dt))){
         idx = (i - 1)*nrow(count_dt) + j
         adj_p_df[i,j] = adj_p[idx]
         adj_s_df[i,j] = assing_significance_label(adj_p[idx])
    }
}

#write.table(p_dt, args[2], row.names = F, col.names = F, quote = F)
#write.table(s_dt, args[3], row.names = F, col.names = F, quote = F)

write.table(adj_p_df, args[2], row.names = F, col.names = F, quote = F)
write.table(adj_s_df, args[3], row.names = F, col.names = F, quote = F)


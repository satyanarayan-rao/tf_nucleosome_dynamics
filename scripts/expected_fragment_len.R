library(pracma)

expected_flen = function(fname){
    dt = read.table(fname, sep="\t", header = F, stringsAsFactors = F)
    
    #print(head(dt))
    cl_ids = unique(dt$V2)
    sum_of_cls = lapply(cl_ids, function(x) {
                            sum(dt[which(dt$V2 == x), "V3"])
    })
    
    names(sum_of_cls) = cl_ids
    prob = c()
    prob_mult_x = c()
    for (i in seq(nrow(dt))){
        #print(i)
        tmp = dt[i, "V3"]/sum_of_cls[[dt[i,"V2"]]]
        #print(tmp)
        prob = c(prob, tmp)
        prob_mult_x = c(prob_mult_x, tmp*dt[i, "V1"])
    }
    dt$prob = prob
    dt$prob_mult_x = prob_mult_x
    #print (head(dt))
    
    expt = lapply(cl_ids, function(x) {
                            sum(dt[which(dt$V2 == x), "prob_mult_x"])
    })
    names(expt) = cl_ids 
    expt_df = data.frame(actual = cl_ids, 
                          expt = unlist(lapply(cl_ids, function (x){
                                         expt[[x]]})))
    expt_df = expt_df[order(expt_df$expt),]
    expt_df$assigned = seq (length(cl_ids))
    #print ("returning")
    return(list(e_df = expt_df, new_dt = dt))
}

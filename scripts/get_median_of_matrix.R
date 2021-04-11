args = commandArgs(trailingOnly = T)
dt = read.table(args[1], sep = "\t", header = F, stringsAsFactors = F)
print (ncol(dt))
head(dt[1:10, 1:10])
dt_sub = dt [1:50, paste0("V", seq(ncol(dt)/2 - 50, ncol(dt)/2 + 50))] 
med = median(as.matrix(dt_sub))
mea = mean (as.matrix(dt_sub))
cat(c(med, mea))

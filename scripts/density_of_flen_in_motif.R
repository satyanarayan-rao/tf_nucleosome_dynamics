# $1: input file : flen count_in_motif
# $2: number of points where density should be calculated
# $3: min_flen: from where density calculation should start
# $4: max_flen: to where density calculation should stop 
# $5: output file
# $6: minimum number of reads
# $7: binwidth
# $8: genomewide background length distribution file: example: metadata/dens_flen_IH02.tsv
args = commandArgs(trailingOnly = T)

input_df = read.table(args[1], sep = "", header = F, stringsAsFactors = F)
density_n = as.integer(args[2])
from = as.integer(args[3])
to = as.integer(args[4])
min_reads = as.integer (args[6])
dt_sub_global = input_df[which(input_df[,6]>=min_reads), seq(5)]
dt_sub = dt_sub_global[, 5, drop = F]
convert = function (vec){as.numeric(unlist(strsplit(vec, split = "-")))}

flen_vec = apply(dt_sub, 1, convert)
names (flen_vec) = row.names (dt_sub)

genomic_background = read.table(args[8], header = F, sep = "\t", stringsAsFactors = F) 
norm_background = genomic_background$V2/sum(genomic_background$V2)
genomic_background$V2 = genomic_background$V2/sum(genomic_background$V2)

expt_genome_len = sum(genomic_background$V1 * genomic_background$V2) 

density_list = lapply(names(flen_vec), function (x) {
                den = density (flen_vec[[x]], bw = as.integer(args[7]),
                               from = from, to = to, n = density_n) 
                #norm_den = den$y/sum(den$y)
                
                #residual_den = norm_den - norm_background
                #div_den = norm_den/norm_background
                return (den$y)}) 
                #return (residual_den)}) 
                #return (div_den)}) 
names(density_list) = names (flen_vec)
density_df = t(round(as.data.frame(density_list), 4))
density_df = cbind (dt_sub_global[, seq(4)], density_df)
write.table(density_df, file = args[5], sep ="\t", row.names = F, col.names = F, quote = F) 

suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("IRanges"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("yaml"))
suppressPackageStartupMessages(library("stringr")) 
suppressPackageStartupMessages(library("R.utils")) 
options(show.error.locations = T)
#options (warn = 2)
source("scripts/utils.R")

args = commandArgs (trailingOnly = T)
tfbs_bed = args[1]
cfdna_bed = args[2]
out_file = args[3]
chunk_size = as.integer (args[4])
len_column = as.integer(args[5])
verbose_mode = as.integer(args[6])
tfbs.data = read.csv(tfbs_bed, 
                     sep = "\t", header = FALSE, 
                     stringsAsFactors=FALSE, comment.char="#")
original_tfbs_order = do.call(paste, c(tfbs.data[, seq(4)], sep = "^"))
names (tfbs.data) = bed_file_header 
tfbs.gr = GRanges (IRanges(start = tfbs.data$`chromStart`,
                           end = tfbs.data$`chromStop`),
                   seqnames = tfbs.data$`chrom`)
tfbs.git = GenomicRanges::GNCList(tfbs.gr) # `git`: genomic interval tree 
cfdna_con = NULL 
if (is_gz(cfdna_bed)){
    #cfdna_con = gzcon(file(cfdna_bed, "rb"))
    cfdna_con = file(cfdna_bed, "rt")
}else{
    cfdna_con = file(cfdna_bed, "r")
}
while_loop_counter = 0
cfdna_file_line_count = as.double (0)
valid_fragments_count = 0 
#print (head (tfbs.gr))
#define a list that will be populated with length of fragment center counts
# for instance, tf_dict[["chr1^1^2^name"]] = "35-76-250"
tf_dict = list ()
#cat (str(tf_dict))
#bt = data.frame(tf_dict)
#print (bt[1:2, 1:50])
#stop()
tf_dict_df = NULL
ver_con = NULL 
if (verbose_mode == 1){
    ver_con = file(paste0(out_file, "verbos.tsv"), "w") 
}
while (TRUE){
    while_loop_counter = while_loop_counter + 1
    #cat (paste0(while_loop_counter, "\n"))
    cfdna.chunk = read_from_connection(con = cfdna_con,
                      chunk_size = chunk_size)
    #print (head(cfdna.chunk))
    if (!is.null(cfdna.chunk)) {
        names(cfdna.chunk) =  bed_file_header
        cfdna_file_line_count = cfdna_file_line_count + nrow(cfdna.chunk)
        #cat (paste0(cfdna_file_line_count, "\n"))
        overlap.idx = find_reads_center_in_tfbs (df = cfdna.chunk,
                          it = tfbs.git)
        cfdna.chunk = cfdna.chunk[overlap.idx@from, ]
        tfbs.it.hits = tfbs.data[overlap.idx@to,]
        #print (head(cfdna.chunk))
        if (nrow(cfdna.chunk) !=0) { # if the hits are found 
            hit_keys = do.call(paste, c(tfbs.it.hits[, seq(4)], sep = "^"))
            if (verbose_mode == 1){
                intersect_df = cbind(tfbs.it.hits, cfdna.chunk)
                write.table(intersect_df,  ver_con, col.names = F, row.names =F, 
                            quote = F, sep = "\t", append = T)
            }
            cnt = 1
            for (k in hit_keys){
                #flen = as.character(cfdna.chunk[cnt, "score"])
                #flen_str = paste(cfdna.chunk[, len_column], collapse="-")
                flen_str = cfdna.chunk[cnt, len_column]
                #tf_dict[[k]][[flen]] = tf_dict[[k]][[flen]] + 1
                if (is.null(tf_dict[[k]])){
                    tf_dict[[k]] = flen_str 
                }else{
                    tf_dict[[k]] = paste(tf_dict[[k]], flen_str, sep="-")
                } 
                cnt = cnt + 1 
            }
           
        } else {
            next
        }

    }else{
        cat ("Finished reading from file", 
             c(cfdna_bed), 
             "closing connection!",  "\n")
        #close(cfdna_con)
        break
    }
}
tf_dict_df = t(as.data.frame(tf_dict))
row.names(tf_dict_df) = names(tf_dict)
rnames = row.names (tf_dict_df)
rnames = str_replace_all(rnames, "\\^", "\t")
tf_dict_df = cbind (tfbs = rnames, tf_dict_df)
no_of_fragments = str_count(tf_dict_df[, 2], "-") 
no_of_fragments = no_of_fragments + 1 
tf_dict_df = cbind(tf_dict_df, flen_count = no_of_fragments)
write.table (tf_dict_df, file = out_file,  
             sep = "\t", col.names = F, row.names = F, quote = F)
if (verbose_mode == 1) {
    close(ver_con)
    gzip (paste0(out_file, "verbos.tsv"), overwrite = T)
}
close (cfdna_con)

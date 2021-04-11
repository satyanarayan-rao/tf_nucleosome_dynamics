suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("IRanges"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("yaml"))
suppressPackageStartupMessages(library("stringr")) 
options(show.error.locations = T)
#options (warn = 2)
source("scripts/utils.R")

args = commandArgs (trailingOnly = T)
tfbs_bed = args[1]
cfdna_bed = args[2]
tfbs.data = read.csv(tfbs_bed, 
                     sep = "\t", header = FALSE, 
                     stringsAsFactors=FALSE, comment.char="#")
names (tfbs.data) = bed_file_header 
tfbs.gr = GRanges (IRanges(start = tfbs.data$`chromStart`,
                           end = tfbs.data$`chromStop`),
                   seqnames = tfbs.data$`chrom`)
tfbs.git = GenomicRanges::GNCList(tfbs.gr) # `git`: genomic interval tree 
cfdna_con = NULL 
if (is_gz(cfdna_bed)){
    cfdna_con = gzcon(file(cfdna_bed, "rb"))
}else{
    cfdna_con = file(cfdna_bed, "r")
}
while_loop_counter = 0
cfdna_file_line_count = as.double (0)
valid_fragments_count = 0 
#print (head (tfbs.gr))
# define a hash with the following structure: 
# key: first four fields of bed entry 
# key->subkeys: 0-700
seed_count = rep(0, 700)
label = as.character(seq(1,700)) 
tf_dict = create_tfbs_hash(df = tfbs.data, 
                                 seed_count = seed_count, label = label)
#cat (str(tf_dict))
#bt = data.frame(tf_dict)
#print (bt[1:2, 1:50])
#stop()
tf_dict_df = NULL
while (TRUE){
    while_loop_counter = while_loop_counter + 1
    #cat (paste0(while_loop_counter, "\n"))
    cfdna.chunk = read_from_connection(con = cfdna_con,
                      chunk_size = 1)
    #print (head(cfdna.chunk))
    if (!is.null(cfdna.chunk)) {
        names(cfdna.chunk) =  bed_file_header
        cfdna_file_line_count = cfdna_file_line_count + nrow(cfdna.chunk)
        overlap.idx = find_reads_center_in_tfbs (df = cfdna.chunk,
                          it = tfbs.git)
        cfdna.chunk = cfdna.chunk[overlap.idx@from, ]
        tfbs.it.hits = tfbs.data[overlap.idx@to,]
        if (nrow(cfdna.chunk) !=0) { # if the hits are found 
            #print (head(overlap.idx))
            #cat ("### --------------------- ### \n")
            #print (head (tfbs.it.hits))
            #cat ("### --------------------- ### \n")
            #print (head (cfdna.chunk)) 
            #break
            print (nrow(cfdna.chunk))
            hit_keys = do.call(paste, c(tfbs.it.hits[, seq(4)], sep = "^"))
            cnt = 1 
            for (k in hit_keys){
                flen = as.character(cfdna.chunk[cnt, "score"])
                tf_dict[[k]][[flen]] = tf_dict[[k]][[flen]] + 1
                cnt = cnt + 1 
            }
            #if (is.null(tf_dict_df)){
            #    tf_dict_df = do.call(rbind, lapply(tf_dict, data.frame)) 
            #}else{
            #    tmp = do.call(rbind, lapply(tf_dict, data.frame))
            #    
            #    tf_dict_df = rbind (tf_dict_df, tmp)
            #}            
           #row.names(tf_dict_df) = names(tf_dict)
            
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
tf_dict_df = do.call(rbind, lapply(tf_dict, data.frame)) 
rnames = row.names (tf_dict_df)
rnames = str_replace_all(rnames, "\\^", "\t")
tf_dict_df = cbind (tfbs = rnames, tf_dict_df)
write.table (tf_dict_df, file = args[3], 
             sep = "\t", col.names = F, row.names = F, quote = F)
close (cfdna_con)

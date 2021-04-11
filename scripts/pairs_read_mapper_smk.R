#  /cluster/software/modules-sw/R/3.5.0/bin/ Rscript
#library("rtracklayer")

suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("IRanges"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("yaml"))
suppressPackageStartupMessages(library("stringr"))

source("scripts/utils.R") # all files have to be relative to the `Snakemake` file

#print (snakemake@params["mapping_params"])
#print (snakemake@params["chunk_size"])
#q()

########### Read inputs from yaml list object coming from snakemake ########### 
params = snakemake@params
mapping_params = params$mapping_params

tfbs.data = read.csv(mapping_params$roi$file, 
                     sep = "\t", header = FALSE, 
                     stringsAsFactors=FALSE, comment.char="#")
tfbs.upstream = 0 - as.integer (mapping_params$roi$left_flank) - 
                    as.integer (mapping_params$roi$length/2 + 0.5) 
names(tfbs.data) = bed_file_header
print(names(tfbs.data))
output.read_mapped.file = snakemake@output$count_matrix
output.read_counts.file = get_other_file_name(output.read_mapped.file, prefix="count")
output.time_log.file = get_other_file_name(output.read_mapped.file, prefix="runtime")

pairs.file = snakemake@input$pairs_list
chunk_size = params$chunk_size
########### Core codes ########### 


# check if output file already exist, if so then delete
if (file.exists(output.read_mapped.file)){
    file.remove(output.read_mapped.file)
}

# create a dataframe with `tfbs.data`'s `chrom`_`chromStart`_`chromStop` as key


tfbs.gr = GRanges (IRanges(start = tfbs.data$`chromStart`,
                           end = tfbs.data$`chromStop`),
                   seqnames = tfbs.data$`chrom`)
                   #strand=tfbs.data$`strand`)
tfbs.git = GenomicRanges::GNCList(tfbs.gr) # `git`: genomic interval tree 

pairs.file.df = read.csv(pairs.file, header = FALSE, sep="", 
                         stringsAsFactors=FALSE, comment.char = "#") 
names(pairs.file.df) = "pairs_files"

valid_fragments_count = 0 
valid_tfbs_region_fragments_count = 0 
pairs_file_line_count = as.double (0)
time1 = Sys.time()
for (pairs.file in pairs.file.df$`pairs_files`){
    con = file(pairs.file, "r")
    while_loop_counter = 0 
    while(TRUE){
        while_loop_counter = while_loop_counter + 1
        df.chunk = read_from_connection(con=con, chunk_size = chunk_size)
        if(!is.null(df.chunk)){
           pairs_file_line_count = pairs_file_line_count + nrow(df.chunk)
            # perform operation on data frame ...
            names(df.chunk) = pairs_file_header
            # applying fragment size filter
            df.chunk = df.chunk[df.chunk$`fragment_length` >= mapping_params$reads$min_length & 
                                df.chunk$`fragment_length` <= mapping_params$reads$max_length, ]
            if (nrow(df.chunk)!=0){
                valid_fragments_count = valid_fragments_count + nrow(df.chunk)
                # find chr position of fragment center point
                overlapped.idx = find_reads_center_in_tfbs(
                                   df = df.chunk, 
                                   it = tfbs.git)
                #print (overlapped.idx)
                df.chunk = df.chunk[overlapped.idx@from, ]
                tfbs.it.hits = tfbs.data[overlapped.idx@to, ]
                #print (c(dim(df.chunk), dim(tfbs.it.hits)))
                if (nrow(df.chunk) !=0){
                    valid_tfbs_region_fragments_count = valid_tfbs_region_fragments_count +
                                                        nrow(df.chunk)
                    # add relative position [-1kb, 1kb] in df.chunk dataframe
                    if (!is.null(mapping_params$`reads`$`fill_in`)){ # have to fill in the reads
                        upstream = mapping_params$`reads`$`fill_in`$`upstream` 
                        downstream = mapping_params$`reads`$`fill_in`$`downstream` 
                        frag_mid_point = as.integer(df.chunk$`fragment_length`/2 + 0.5)
                        mid_point = tfbs.upstream + 
                                    (df.chunk$`chromStart` - tfbs.it.hits$`chromStart`) + 
                                    frag_mid_point 
                        fill_in_str = paste (as.character(mid_point - min(upstream, frag_mid_point - 1 )),
                                             as.character(mid_point + min(downstream, frag_mid_point - 1)),
                                             sep = "|" )
                        #print (dim(df.chunk))
                        #print (dim(tfbs.it.hits))
                        df.chunk["relative_pos"] = fill_in_str  
                                                
                        
                    } else{
                        upstream = 0  
                        downstream = 0  
                        mid_point = tfbs.upstream +
                               (df.chunk$`chromStart` - tfbs.it.hits$`chromStart`) +
                               as.integer(df.chunk$`fragment_length`/2 + 0.5) 
                        fill_in_str = paste (as.character(mid_point - upstream),
                                             as.character(mid_point + downstream),
                                             sep = "|" )
                        df.chunk["relative_pos"] = fill_in_str  
                    }
                        
                    #df.chunk["tfbs_left_1kb_flank"] = tfbs.it.hits@ranges@start
                    df.chunk["tfbs_left_flank"] = tfbs.it.hits$`chromStart`
                    df.chunk["tfbs_right_flank"] = tfbs.it.hits$`chromStop`
                    data.table::fwrite(df.chunk, output.read_mapped.file,  
                                       sep="\t", row.names=FALSE, append=TRUE)
                } else{
                    next 
                } 
           } else{
               next
           }
    
        } else{
            cat ("Finished reading from file", 
                 c(pairs.file), 
                 "closing connection!",  "\n")
            close(con)
            break
        }
        cat ("finished processing", 
             c(while_loop_counter * chunk_size), 
             "lines of input file", c(pairs.file), "\n" )
    }
}

time2 = Sys.time()
s = difftime (time2, time1, tz, units = c("mins"))
to_write = capture.output(s)

# write count statistics 

count_statistics = list ("total_valid_fragments" = valid_fragments_count,
                         "total_valid_tfbs_region_framgments" = valid_tfbs_region_fragments_count,
                         "pairs_file_line_count" = pairs_file_line_count)

count_statistics = data.frame(count_statistics)
data.table::fwrite(count_statistics,
                   output.read_counts.file, 
                   sep="\t", row.names=FALSE)

# write time taken to run the program

write.csv(to_write, quote=FALSE, row.names=FALSE, file=output.time_log.file)

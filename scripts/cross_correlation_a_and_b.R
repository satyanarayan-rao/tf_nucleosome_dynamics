library(data.table)
library(readr)
library(dplyr)
library(R.utils)
library(tools)
options(error=traceback)
args = commandArgs(trailingOnly=TRUE)

file_x = args[1]
file_y = args[2]
lag.val = as.numeric(args[3])
out_file_name = tools::file_path_sans_ext(args[4])
nan_entries_file = args[5]
out_file_con = file(out_file_name, "w")
out_file_nan_con = file(nan_entries_file, "w")

con_x = file (file_x, "r")
con_y = file (file_y, "r")

while (TRUE){
    line_x = readLines(con_x, n = 1)
    line_y = readLines(con_y, n = 1)
    if (length (line_x) == 0 || length(line_y) == 0){
        break
    }
    # parse line data and make a float vector
    x_t_with_name = unlist (strsplit(line_x, "\t"))
    x_t = as.numeric (x_t_with_name [2:length(x_t_with_name)]) 
    row_name = x_t_with_name[1]

    y_t_with_name = unlist (strsplit(line_y, "\t"))
    y_t = as.numeric (y_t_with_name [2:length(y_t_with_name)]) 
    # do cross-correlation 
    x_y_ccf = ccf (x_t, y_t, plot =F, lag.max = lag.val, na.action = na.pass) 
    if (all(is.na(x_y_ccf$acf))) {
      nan_entry_str = row_name
      writeLines(nan_entry_str, out_file_nan_con)
    } else {
      #print (x_y_ccf$acf)
      ccf_s = paste(round(x_y_ccf$acf, 3), collapse="\t")
      # store the result in a string vector
      ccf_str = paste(row_name, ccf_s, collpase = "\t")
      #print (ccf_str)
      #ccf_str = paste (ccf_str, "\n")
      writeLines(ccf_str, out_file_con)
    }
}
close (out_file_con)
close (out_file_nan_con)
gzip(out_file_name, destname=paste0(out_file_name, ".gz"), overwrite=T) 

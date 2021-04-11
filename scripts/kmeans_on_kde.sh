#!/bin/bash
# $1: input kde file
# $2: nclust
# $3: max_iter
# $4: output file
Rscript scripts/kmeans_of_kde_matrix.R $1 $2 $3 $4 $5


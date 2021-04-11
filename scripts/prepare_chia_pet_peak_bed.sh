#!/bin/bash
# $1: peak_location_tsv
# $2: peak_location_bed

# $1 has tab delimiter

awk -F'[-:\t]' '{OFS="\t"} {print $1, $2, $3, $1":"$2"-"$3"|"$4"|"$5"|"$NF}'  $1 > $2

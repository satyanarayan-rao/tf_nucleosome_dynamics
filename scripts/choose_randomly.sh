#!/bin/bash
# $1: observed
# $2: whole random set
# $3: output

total_to_extract=`wc -l $1 | awk '{print $1}'`
cat $2 | shuf -n ${total_to_extract} > $3

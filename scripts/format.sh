#!/bin/bash
awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"|"lab"`"NR}' lab=$3 $1 > $2 

from __future__ import division
import sys
import gzip
import getopt
import re
import math
import numpy as np
import scipy
from scipy.signal import find_peaks
from numpy import pi
from scipy.signal import savgol_filter
from collections import defaultdict
import math

peak_dict = {} 
input_fp = gzip.open(sys.argv[1])
slop_value = int(sys.argv[2])
out_fp = gzip.open(sys.argv[3], "wb")
bed_with_dnase_loc_fp = open(sys.argv[4], "w")
min_of_all = float('Inf')
max_of_all = 0 - float('Inf')
min_of_peaks = float('Inf')
max_of_peaks = 0 - float ('Inf')
eom_dict = {}
site_to_order = []
peaks_to_order = []
nearest_peak_with_sign  = {}
in_or_out_dict = {}   
#cols_in_file = len(bytes.decode(input_fp.readline()).split("\t")) - 1 # first  col is chr name
#print(cols_in_file)
#mid_point = cols_in_file/2 

#input_fp.seek(0,0)
for line in input_fp: 
    line_txt = bytes.decode(line)
    first_tab = line_txt.find("\t")
    site_id = line_txt[0:first_tab]
    signal_vec = np.fromstring(line_txt[first_tab + 1:], sep = "\t", dtype = float)
    mid_point = len (signal_vec)/2  
    peak_left_boundary = slop_value 
    peak_right_boundary = len(signal_vec) - slop_value
    width = len(signal_vec) 
    eom_dict[site_id] = list(signal_vec)
    smooth_signal = savgol_filter(signal_vec, 9, 1)
    peaks, _ = find_peaks(smooth_signal, 
             height=3*smooth_signal.std()+smooth_signal.mean(), distance=25)
    peak_dict[site_id] = peaks
    site_to_order.append(site_id)
   
    if len(peaks) == 0: 
        peaks_to_order.append(float('Inf'))
        nearest_peak_with_sign[site_id] = float('Inf')
        in_or_out_dict[site_id] = 1 
    else: 
        min_peak = min(peaks)
        peaks_to_order.append(abs(min(peaks) - mid_point))
        nearest_peak_with_sign[site_id] = min_peak - mid_point
        if (min_peak < peak_left_boundary) or (min_peak > peak_right_boundary):
            in_or_out_dict[site_id] = 1 
        else:
            in_or_out_dict[site_id] = 0 
        #print (nearest_peak_with_sign[site_id])
        min_of_peaks = min(peaks) 
        max_of_peaks = max(peaks)
    if min_of_all > min_of_peaks: 
        min_of_all = min_of_peaks
    if max_of_all < max_of_peaks:
        max_of_all = max_of_peaks

print([min_of_all,max_of_all])
idx_sorted_peaks = sorted(range(len(peaks_to_order)), 
                     key = lambda k : peaks_to_order[k])[::-1] 
sorted_sites = [site_to_order[i] for i in idx_sorted_peaks]
#min_with_offset = max(0, min_of_all - 25)
#max_with_offset = min(2000, min_of_all + 25 )
for s in sorted_sites:
    to_write = s + "\t" + "\t".join(map(str, eom_dict[s])) 
    bytes_to_write = bytes(to_write + "\n", encoding = "ascii")
    out_fp.write(bytes_to_write) 
    bed_str = "\t".join(s.split("@"))
    nearest_peak = str(nearest_peak_with_sign[s])  
    in_or_out = str(in_or_out_dict[s])
    to_write_bed = bed_str + "\t" + nearest_peak + "\t" + in_or_out
    bed_with_dnase_loc_fp.write(to_write_bed + "\n")

   
input_fp.close()
bed_with_dnase_loc_fp.close()

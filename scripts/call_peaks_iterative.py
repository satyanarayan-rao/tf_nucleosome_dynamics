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


def readwig (inFile):

##bedGraph section chr1:0-828999
#chr1    0    9999    0
#chr1    9999    10099    12.516
#chr1    10099    10199    15.0193
#chr1    10199    10299    1.2516
    wig=open(inFile, 'r' )
    val={}
    val_ar=[]
    for line in wig:
        lineL = line.split()
        if "bedGraph" in line:
            cj=lineL[-1]
            chrom=cj.split(":")[0].replace("chr", "")
            #print(chrom)
            val.setdefault(chrom,{})
        elif "rack" not in line:
            st = int (lineL[1])
            en = int (lineL[2])
            pos = int((st + en)/2)
            val[chrom][pos] = float(lineL[-1])
            val_ar.append(float(lineL[-1]))
    return(val,val_ar)    



def main(argv):
    
    opts, args = getopt.getopt(argv[1:], "i:o:c:", ["infile=", "outfile=","cutoff="])
    
    if(len(argv)<2):
        print("\nusage: python %s -i/--infile -o/--outfile\n" % argv[0])
        sys.exit(2)
    
    infile=""
    outfile=""
    cutoff=0
    for opt, arg in opts:
        print(opt + " " + arg)
        if opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-o", "--outfile"):
            outfile = arg
        elif opt in ("-c", "--cutoff"):
            cutoff = int(arg)
    profile, val_ar = readwig(infile)
    nval_ar = np.array(val_ar)
    
    #Get genome-wide cutoff iteratively
    
    smooth_ar = savgol_filter(np.array(nval_ar),5,1)
    smooth_ar = np.sort ( smooth_ar )[::-1]
    
    npeaks=0
    prev_peaks=-1
    it=0
    it_cutoff=0

    while (npeaks>prev_peaks):
        print(str(it)+" Npeak: "+str(npeaks)+" Prev_peak: "+str(prev_peaks)+" N: "+str(smooth_ar.size)+" Cutoff: "+str(it_cutoff))
        prev_peaks = npeaks
        cont = []
        m = smooth_ar.mean()
        s = smooth_ar.std()
        for i in np.arange(smooth_ar.size):
            val = (smooth_ar[i]-m)/s
            #print(val)
            if(val >= cutoff):
                it_cutoff = smooth_ar[i]
                cont.append(i)
                npeaks+=1
            else:
                break
            smooth_ar = np.delete(smooth_ar,cont)
            it+=1
    
    print(str(it)+" Npeak: "+str(npeaks)+" Prev_peak: "+str(prev_peaks)+" N: "+str(smooth_ar.size)+" Cutoff: "+str(it_cutoff))

    #For each chromosome, call peaks
    fh = open(outfile,"w")
    for j in profile.keys(): #Every chromosome
        temp_list = [] #Values
        pos_list = [] #Chromosome positions
        for i in sorted(profile[j].keys()): #Every position
            temp_list.append(profile[j][i])
            pos_list.append(i)
        chrom_array = np.array(temp_list)
        #print(chrom_array.ndim)
        #print(len(chrom_array))
        smooth_ar = savgol_filter(chrom_array,5,1)
        peaks, _ = find_peaks(smooth_ar, height=it_cutoff, distance=2)
        
        for p in peaks:
            peak_pos = pos_list[p]
            peak_val = temp_list[p]
            temp_str = j + "\t" + str(peak_pos) + "\t" + str(peak_pos+1) + "\t" + j + "_" + str(peak_pos) + "\t" + str(peak_val) + "\n"
            fh.write(temp_str)
    
    fh.close()
if __name__ == "__main__":
    main(sys.argv)

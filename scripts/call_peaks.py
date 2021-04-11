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
	wig=open(inFile, 'r' )
	val={}
	val_ar=[]
	for line in wig:
		lineL = line.split()
		if "chr" in lineL[1]:
			cj=lineL[1]
			chrom=cj[9:]
			print(chrom)
			val.setdefault(chrom,{})
		elif "rack" not in line:
			pos=int(lineL[0])
			val[chrom][pos] = float(lineL[1])
			val_ar.append(float(lineL[1]))
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
			print("In Here")
			infile = arg
		elif opt in ("-o", "--outfile"):
			print("Out Here")
			outfile = arg
			print(outfile)
		elif opt in ("-c", "--cutoff"):
			cutoff = int(arg)
	profile, val_ar = readwig(infile)
	nval_ar = np.array(val_ar)
	#For each chromosome, call peaks
	fh = open(outfile,"w")
	for j in profile.keys(): #Every chromosome
		temp_list = [] #Values
		pos_list = [] #Chromosome positions
		for i in sorted(profile[j].keys()): #Every position
			temp_list.append(profile[j][i])
			pos_list.append(i)
		chrom_array = np.array(temp_list)
		print(chrom_array.ndim)
		print(len(chrom_array))
		smooth_ar = savgol_filter(chrom_array,9,1)
		peaks, _ = find_peaks(smooth_ar, height=cutoff*nval_ar.std()+nval_ar.mean(), distance=25)
		#for i in np.arange(0,len(smooth_ar)):
		#	print(str(pos_list[i])+"\t"+str(smooth_ar[i]))
		
		for p in peaks:
			peak_pos = pos_list[p]
			temp_str = j + "\t" + str(peak_pos-250) + "\t" + str(peak_pos+250) + "\t" + str(p) + "\n"
			fh.write(temp_str)
	
	fh.close()
if __name__ == "__main__":
    main(sys.argv)

import os
import re
import glob
import numpy as np
from scipy.stats import scoreatpercentile                                                                         
from numpy import sort      

def concatenate1(regex):
	for f in glob.glob(regex):
		os.system("cat "+f+" >> out1.dat")

def concatenate2(regex):
	for f in glob.glob(regex):
		os.system("cat "+f+" >> out2.dat")

def convert_to_csv(inp, outp):
	with open(inp, 'r') as infile:
	    with open(outp, 'w') as outfile:
	        regex =re.compile(r'#')
	        for line in infile:
	            columns = line.strip().split()
	            if regex.match(line) is None:
	                outfile.write(",".join(columns)+"\n")

	            # if regex.match(line) is None or np.isnan(line) is False or np.isinf(line) is False:


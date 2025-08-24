#!/usr/bin/env python

import os,sys
from string import *

import scipy, scipy.stats

infile = sys.argv[1]
outfile = sys.argv[2]
f = open(infile, 'r')
outf = open(outfile, 'w')
outf.write(f.readline())
for line in f:
	tline = line.strip()
	sline = tline.split('\t')
	gn = sline[0]
	cv = []
	for i in range(len(sline)-1):
		cv.append(float(sline[i+1]))
	
	mtcv = scipy.stats.zscore(cv)
	outf.write(gn)
	for c in mtcv:
		outf.write("\t%f" % c)
	outf.write('\n')
outf.close()
f.close()

		
		
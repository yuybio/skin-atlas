#!/usr/bin/env python

import os,sys
from string import *

import scipy, scipy.stats
import statistics

f = open(sys.argv[1], 'r')
headerl = f.readline()
tline = headerl.strip()
sline = tline.split('\t')
samples=[]
sample_expv={}
for i in range(len(sline)-1):
    samples.append(sline[i+1])
    v = []
    sample_expv.update({sline[i+1]:v})



for line in f:
    tline=line.strip()
    sline = tline.split('\t')
    gn = sline[0]
    for i in range(len(sline)-1):
        s = samples[i]
        sample_expv[s].append(float(sline[i+1]))
f.close()

outf=open(sys.argv[2], 'w')
outf.write("sample\tsite\tz\n")
for s in sample_expv.keys():
    m = statistics.mean(sample_expv[s])
    site = s[0:2]
    outf.write("%s\t%s\t%f\n" %(s, site, m))
outf.close()

        
#!/usr/bin/env python

import os,sys
from string import *
import scipy.stats

infile=sys.argv[1]
outfile=sys.argv[2]

f = open("module_exp_matrix.txt", 'r')

exp = {}
header = f.readline()
for line in f:
    tline = line.strip()
    sline = tline.split('\t')
    if (len(sline) != 71):
        print("error column size %d\t%s\n" % (len(sline), sline[0]) )
        exit(1)
    ct = sline[0]
    
    exp.update({ct:line})
f.close()

f = open("celltype_sample_proportion.txt", 'r')

f.readline()
for line in f:
    tline = line.strip()
    sline = tline.split('\t')
    if (len(sline) != 71):
        print("error column size %d\t%s\n" % (len(sline), sline[0]) )
        exit(1)
    ct = sline[0]
    cv = []
    exp.update({ct:line})
f.close()
f = open("granular_celltype_sample_proportion.txt", 'r')
f.readline()
for line in f:
    tline = line.strip()
    sline = tline.split('\t')
    if (len(sline) != 71):
        print("error column size %d\t%s\n" % (len(sline), sline[0]) )
        exit(1)
    ct = sline[0]
    exp.update({ct:line})
f.close()

f = open(infile, 'r')
f.readline()
terms = []
for line in f:
    tline = line.strip()
    sline = tline.split('\t')
    term = sline[0]
#    if term not in exp.keys:
 #       print("error cannot find term %s \n"%term )
 #       exit(1)
    terms.append(term)
f.close()

outf = open(outfile, 'w')
outf.write(header)
for t in terms:
    outf.write(exp[t])
outf.close()

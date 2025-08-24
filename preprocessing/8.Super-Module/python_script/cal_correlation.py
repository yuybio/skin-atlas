#!/usr/bin/env python

import os,sys
from string import *
import scipy.stats

f = open("module_exp_matrix.txt", 'r')

exp = {}
f.readline()
for line in f:
    tline = line.strip()
    sline = tline.split('\t')
    if (len(sline) != 71):
        print("error column size %d\t%s\n" % (len(sline), sline[0]) )
        exit(1)
    ct = sline[0]
    cv = []
    for i in range(1,len(sline)):
        c = float(sline[i])
        cv.append(c)
    print("%d\n" %len(cv) )
    exp.update({ct:cv})
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
    for i in range(1,len(sline)):
        c = float(sline[i])
        cv.append(c)
    print("%d\n" %len(cv) )
    exp.update({ct:cv})
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
    cv = []
    for i in range(1,len(sline)):
        c = float(sline[i])
        cv.append(c)
    print("%d\n" %len(cv) )
    exp.update({ct:cv})
f.close()

def remove_zero( ve1, ve2 ):
    rve1 = []
    rve2 = []
    for i in range(0,len(ve1)):
        if ( (ve1[i] !=0) and (ve2[i] !=0) ):
            rve1.append(ve1[i])
            rve2.append(ve2[i])
    return (rve1, rve2)

outf1 = open("moduleexp_composition_correlation.txt", 'w')
outf2 = open("moduleexp_composition_correlation_p.txt", 'w')

for ct in exp.keys():
    outf1.write("\t%s" % ct)
    outf2.write("\t%s" % ct)
outf1.write("\n")
outf2.write("\n")

for ct1 in exp.keys():
    outf1.write("%s" %ct1)
    outf2.write("%s" %ct1)
    
    cv1 = exp[ct1]
    for ct2 in exp.keys():

        cv2 = exp[ct2]
        (rve1, rve2) = remove_zero(cv1, cv2)
        r, p = scipy.stats.pearsonr(rve1, rve2)
        outf1.write("\t%f" %r)
        outf2.write("\t%f" %p)
    outf1.write("\n")
    outf2.write("\n")
outf1.close()
outf2.close()


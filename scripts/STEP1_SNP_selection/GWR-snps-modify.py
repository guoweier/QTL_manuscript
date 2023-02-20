#! /usr/bin/env python

import os, sys, math, time
from optparse import OptionParser

usage = "\n\n%prog"
parser = OptionParser(usage=usage)

parser.add_option("-s", "--snp", dest="s", help="Input SNP File")
parser.add_option("-o", "--outputfile", dest="o", help="Output file.")

(opt, args) = parser.parse_args()

snpfile=open(opt.s,"r")
o=open(opt.o,"w")

SNPhead = snpfile.readline()
modify = []

for line in snpfile:
    ls=line.split('\t')
    male=ls[2].split('_')[1]
    female=ls[4].split('_')[1]
    if len(male) == 1:
        if len(female) == 1:
            modify.append(ls)
        else:
            fhet=female.split('/')
            if fhet[0] != '*' and fhet[1] != '*':
                modify.append(ls)
    else:
        mhet=male.split('/')
        if mhet[0] != '*' and mhet[1] != '*':
            if len(female) == 1:
                modify.append(ls)
            else:
                fhet=female.split('/')
                if fhet[0] != '*' and fhet[1] != '*':
                    modify.append(ls)

snpfile.close()

o.write(SNPhead)
for each in modify:
    ln=each[0]+'\t'+each[1]+'\t'+each[2]+'\t'+each[3]+'\t'+each[4]+'\t'+each[5]
    o.write(ln)

o.close()

            

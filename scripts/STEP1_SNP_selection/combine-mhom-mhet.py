#! /usr/bin/env python

import os, sys, math, time
from optparse import OptionParser

usage = "\n\n%prog"
parser = OptionParser(usage=usage)

parser.add_option("-a", "--mhom", dest="mhom", help="Input male homozygous SNP File")
parser.add_option("-b", "--mhet", dest="mhet", help="Input male heterozygous SNP File")
parser.add_option("-o", "--outputfile", dest="o", help="Output file.")

(opt, args) = parser.parse_args()

homfile=open(opt.mhom, "r")
hetfile=open(opt.mhet, "r")
o=open(opt.o, "w")

hompos=[]
hetpos=[]

for line in homfile:
    ls=line.split('\t')
    if ls[0][0] == 'C':
        hompos.append(ls)

for line in hetfile:
    ls=line.split('\t')
    if ls[0][0] == 'C':
        hetpos.append(ls)

homfile.close()
hetfile.close()

hompos.append(["Chr",float("inf")])
hetpos.append(["Chr",float("inf")])
i=0
j=0
combine=[]
for k in range(len(hompos)-1+len(hetpos)-1):
    if float(hompos[i][1]) <= float(hetpos[j][1]):
        combine.append(hompos[i])
        i=i+1
    elif float(hompos[i][1]) > float(hetpos[j][1]):
        combine.append(hetpos[j])
        j=j+1
header="Chrom\tPos\tMale\tcov-Male\tFemale\tcov-Female"
o.write(header+'\n')

for each in combine:
    ln = each[0]+'\t'+each[1]+'\t'+each[2]+'\t'+each[3]+'\t'+each[4]+'\t'+each[5]
    o.write(ln)

o.close()


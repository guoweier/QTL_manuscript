import os, sys, math, time
from optparse import OptionParser

usage = "\n\n%prog"
parser = OptionParser(usage=usage)

parser.add_option("-H", "--Haplotype", dest="H", help="Input Haplotype File")
parser.add_option("-m", "--outputmale", dest="m", help="Output male file")
parser.add_option("-f", "--outputfemale", dest="f", help="Output female file")

(opt, args) = parser.parse_args()

#haplotype file
haplofile = open(opt.H)
#male and female output file
om = open(opt.m, 'w')
of = open(opt.f, 'w')

haplohead = haplofile.readline()
headerm = 'Chrom\tPos\tM1\tM2'
headerf = 'Chrom\tPos\tF1\tF2'
om.write(headerm+'\n')
of.write(headerf+'\n')

for line in haplofile:
    ln = line.split('\t')
    chrom = ln[0]
    pos = ln[1]
    m1 = ln[2]
    m2 = ln[3]
    f1 = ln[4]
    f2 = ln[5]
    if m2 == '.':
        text1 = chrom+'\t'+pos+'\t'+m1+'\t'+m1
        om.write(text1+'\n')
    if f2 == '.':
        text2 = chrom+'\t'+pos+'\t'+f1+'\t'+f1
        of.write(text2+'\n')
    
om.close()
of.close()
haplofile.close()

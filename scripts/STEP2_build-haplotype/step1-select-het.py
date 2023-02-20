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
om.write(haplohead)
of.write(haplohead)

for line in haplofile:
    ln = line.split('\t')
    m1 = ln[2]
    m2 = ln[3]
    f1 = ln[4]
    f2 = ln[5]
    if m2 != '.':
        om.write(line)
    if f2 != '.':
        of.write(line)
    
om.close()
of.close()
haplofile.close()

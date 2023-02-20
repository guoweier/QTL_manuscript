import os, sys, math, time
from optparse import OptionParser

usage = "\n\n%prog"
parser = OptionParser(usage=usage)

parser.add_option("-a", "--allele", dest="s", help="Input unfiltered call alleles File")
parser.add_option("-p", "--percent", dest="p", type="float", default=0.9, help="the percent of samples to have coverage over the expected number (default is 90%).") 
parser.add_option("-c", "--coverage", dest="c", type="float", default=20.0, help="Input expected total coverage reads number for each position (default is 20 reads).")
parser.add_option("-o", "--outputfile", dest="o", help="Output file.")

(opt, args) = parser.parse_args()

allelefile = open(opt.s)
expcov = opt.c
per = opt.p
o = open(opt.o,'w')

indices = [] #the position of total coverage column for each sample

allelehead = allelefile.readline()
ah = allelehead.split('\t')
if allelehead[0] == 'C':
    count = 0
    for m in ah:
        if m.split('-')[0] == 'TotalCov':
            indices.append(count)
        count += 1
print(indices)
o.write(allelehead)

samplenum = int(len(indices)) #the number of samples that need to have expected coverage number
print(samplenum)
for line in allelefile:
    satcov = 0
    ln = line.split('\t')
    for i in indices:
        if ln[i] == '.':
            continue
        elif int(ln[i]) >= expcov:
            satcov += 1
    if satcov >= samplenum*per:
        o.write(line)

o.close()
allelefile.close()


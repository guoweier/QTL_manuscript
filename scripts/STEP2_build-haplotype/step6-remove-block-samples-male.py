import os, sys, math, time
from optparse import OptionParser
from collections import OrderedDict

usage = "\n\n%prog"
parser = OptionParser(usage=usage)

parser.add_option("-f", "--file", dest="f", help="input removed weird SNP file")
parser.add_option("-s", "--selectfile", dest="s", help="input file with heterozygous SNPs and RNA-seq samples")
parser.add_option("-o", "--outputfile", dest="o", help="output file")

(opt, args) = parser.parse_args()

cfile = open(opt.f) #open blocked haplotype list
sfile = open(opt.s) 
o = open(opt.o, "w") #output weird SNP removed file

ftb = [] #table for store removed weird SNPs
outtb = [] #table for final output SNPs

chead = cfile.readline()
for line in cfile:
	ln = line.split("\t")
	item = [ln[0],ln[1]]
	ftb.append(item)

shead = sfile.readline()
o.write(shead)

for line in sfile:
	ln = line.split("\t")
	snp = [ln[0],ln[1]]
	if snp in ftb:
		o.write(line)
	else:
		continue

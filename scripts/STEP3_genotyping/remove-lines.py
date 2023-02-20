import os, sys, math, time
from optparse import OptionParser

usage = "\n\n%prog"
parser = OptionParser(usage=usage)

parser.add_option("-s", "--snp", dest="s", help="Input generated haplotype file")
parser.add_option("-o", "--outputfile", dest="o", help="Output file.")

(opt, args) = parser.parse_args()

f = open(opt.s)
o = open(opt.o, "w")

for line in f:
	ln = line.split("\t")
	if ln[0] != "-":
		o.write(line)

f.close()
o.close()
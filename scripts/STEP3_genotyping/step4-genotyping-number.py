import os, sys, math, time
from optparse import OptionParser

usage = "\n\n%prog"
parser = OptionParser(usage=usage)

parser.add_option("-g", "--genotype", dest="g", help="Input step3-genotyping-samples output file")
parser.add_option("-f", "--female", dest="f", type="float", default=1.0, help="whether the input marker list is male or female haplotypes. Male=1, Female=2. Default is male")
parser.add_option("-o", "--outputfile", dest="o", help="Output file.")

(opt, args) = parser.parse_args()

genofile = open(opt.g)
o = open(opt.o, "w")
mf = opt.f

header = genofile.readline()
o.write(header)
for line in genofile:
	line = line[:-1]
	ln = line.split("\t")
	m1 = ln[2]
	m2 = ln[3]
	f1 = ln[4]
	f2 = ln[5]
	text = ln[0]+"\t"+ln[1]+"\t"+m1+"\t"+m2+"\t"+f1+"\t"+f2
	if mf == 1:
		for i in range(6,len(ln)):
			if ln[i] == m1:
				text += "\t"+"1"
			elif ln[i] == m2:
				text += "\t"+"2"
			else:
				text += "\t"+"1.5"
		o.write(text+"\n")
	elif mf == 2:
		for i in range(6,len(ln)):
			if ln[i] == f1:
				text += "\t"+"1"
			elif ln[i] == f2:
				text += "\t"+"2"
			else:
				text += "\t"+"1.5"
		o.write(text+"\n")
genofile.close()
o.close()
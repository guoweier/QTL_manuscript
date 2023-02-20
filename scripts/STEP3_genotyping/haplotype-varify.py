import os, sys, math, time
from optparse import OptionParser
from collections import OrderedDict

#this script is to check for the correctness of generated haplotype list.


usage = "\n\n%prog"
parser = OptionParser(usage=usage)

parser.add_option("-f", "--file", dest="f", help="input file for hapotype with RNA-seq samples results")
parser.add_option("-p", "--percentage", dest="p", type= "float", help="percentage deciding haplotype transfer")
parser.add_option("-n", "--samplenum", dest="n", type="int", help="total sample numbers")
parser.add_option("-o", "--outputfile", dest="o", help="output file")

(opt, args) = parser.parse_args()

hfile = open(opt.f)
o = open(opt.o, "w")
per = opt.p
num = opt.n

header = hfile.readline()
ftable = []
trans = {}
err = []
temp = hfile.tell()

for line in hfile:
	line = line.split("\n")[0]
	ln = line.split("\t")
	ftable.append(ln)


for i in range(len(ftable)):
	ni = i+1
	chrom = int(ftable[i][0][3:])
	if chrom not in trans:
		trans[chrom] = {}
	if ni < len(ftable) and ftable[ni][0] == ftable[i][0]:
		trans[chrom][ftable[i][1]] = 0
		for j in range(6,len(ftable[ni])):
			h1 = ftable[i][j]
			h2 = ftable[ni][j]
			try:
				if float(h1) != float(h2) and float(h1) != 1.5 and float(h2) != 1.5:
					trans[chrom][ftable[i][1]] += 1
			except:
				continue
	elif ni < len(ftable) and ftable[ni][0] != ftable[i][0]:
		continue
otrans = {}
for g in trans:
	otrans[g] = OrderedDict(sorted(trans[g].items()))


for chrom in otrans:
	if len(str(chrom)) == 1:
		chrm = "Chr"+"0"+str(chrom)
	elif len(str(chrom)) == 2:
		chrm = "Chr"+str(chrom)
	for pos in otrans[chrom]:
		if otrans[chrom][pos] >= per*num:
			err.append([chrm, pos])

ohead = "Chrom"+"\t"+"Pos"
o.write(ohead+"\n")
for item in err:
	text = item[0]+"\t"+item[1] 
	o.write(text+"\n")

o.close()


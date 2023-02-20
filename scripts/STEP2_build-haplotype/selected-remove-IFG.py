import os, sys, math, time
from optparse import OptionParser

usage = "\n\n%prog"
parser = OptionParser(usage=usage)

parser.add_option("-s", "--selected", dest="s", help="Input selected het haplotype file")
parser.add_option("-o", "--outputfile", dest="o", help="Output file.")

(opt, args) = parser.parse_args()

#open files#
f = open(opt.s)
o = open(opt.o, "w")

indices = []
count = 0

fhead = f.readline()
fheader = fhead.split("\t")
for item in fheader:
	if item.split("-")[0] in ["Snptype", "SNP1", "SNP2", "%M", "TotalCov", "CovM"]:
		if item.split("-")[1] != "IFG" and item.split("-")[3] != "47":
			indices.append(count)
	count += 1

#output file header#
header = "Chrom\tPos\tM1\tM2\tF1\tF2"
for a in indices:
    factor = fheader[a]
    header += "\t"+factor
o.write(header)

for line in f:
    ln = line.split("\t")
    wline = ln[0]+"\t"+ln[1]+"\t"+ln[2]+"\t"+ln[3]+"\t"+ln[4]+"\t"+ln[5]
    for a in indices:
        factor = ln[a]
        wline += "\t"+factor
    o.write(wline)

f.close()
o.close()

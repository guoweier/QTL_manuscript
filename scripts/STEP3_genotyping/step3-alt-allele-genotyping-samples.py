from __future__ import division
import os, sys, math, time
from optparse import OptionParser
from collections import OrderedDict

# This script uses parent marker haplotype list to genotyping every RNA sequencing sample.
# Different from the previous Step3 script, this one only select SNPs that only inherited alternative alleles. 

usage = "\n\n%prog"
parser = OptionParser(usage=usage)

parser.add_option("-a", "--allele", dest="a", help="Input step2-haplotype-call-allele output file")
parser.add_option("-s", "--sample", dest="s", help="Select the sample name, i.e. GWR-100-255")
parser.add_option("-f", "--female", dest="f", type="float", default=1.0, help="whether the input marker list is male or female haplotypes. Male=1, Female=2. Default is male")
parser.add_option("-o", "--outputfile", dest="o", help="Output file.")

(opt, args) = parser.parse_args()

#open files and parameters
afile = open(opt.a)
sample = opt.s
o = open(opt.o, "w")
mf = opt.f


head = afile.readline()
head = head.split("\n")[0]
hd = head.split("\t")

#write output file header
ohead = "Chrom"+"\t"+"Pos"+"\t"+"M1"+"\t"+"M2"+"\t"+"F1"+"\t"+"F2"+"\t"+sample
o.write(ohead+"\n")



#record the position of selected sample
count = 0
sampletb = []
for item in hd:
	if sample in item and item.split("-")[3] == sample.split("-")[2] and item.split("-")[1] == sample.split("-")[0]:
		if item.split("-")[0] in ["SNP1","SNP2","TotalCov"]:
			sampletb.append(count)
	count += 1

#start read lines from input file
alt_dic = {}

for line in afile:
	line = line.split("\n")[0]
	ln = line.split("\t")
	chrom = ln[0]
	pos = int(ln[1])
	m1 = ln[2]
	m2 = ln[3]
	f1 = ln[4]
	f2 = ln[5]
	snp1 = ln[sampletb[0]]
	snp2 = ln[sampletb[1]]
	cov = ln[sampletb[2]]

	if snp2 == ".":
		continue
	elif snp2 != ".":
		if chrom not in alt_dic:
			alt_dic[chrom] = {}

		#male haplotype
		if mf == 1.0:
			if f2 == ".":
				if m1 == f1:
					alt = m2
					if snp1 == alt or snp2 == alt:
						alt_dic[chrom][pos] = [m1,m2,f1,f2,str(2)]
				elif m2 == f1:
					alt = m1
					if snp1 == alt or snp2 == alt:
						alt_dic[chrom][pos] = [m1,m2,f1,f2,str(1)]
				elif m1 != f1 and m2 != f1:
					alt = m1
					alt_dic[chrom][pos] = [m1,m2,f1,f2,str(1)]
			elif f2 != ".":
				if m1 != f1 and m1 != f2:
					alt = m1
					if snp1 == alt or snp2 == alt:
						alt_dic[chrom][pos] = [m1,m2,f1,f2,str(1)]
				elif m2 != f1 and m2 != f2:
					alt = m2
					if snp1 == alt or snp2 == alt:
						alt_dic[chrom][pos] = [m1,m2,f1,f2,str(2)]
				elif f1 != m1 and f1 != m2 and f2 != m1 and f2 != m2:
					alt = m1
					if snp1 == alt or snp2 == alt:
						alt_dic[chrom][pos] = [m1,m2,f1,f2,str(1)]
			else:
				alt = "."
		
		#female haplotype
		elif mf == 2.0:
			if m2 == ".":
				if f1 == m1:
					alt = f2
					if snp1 == alt or snp2 == alt:
						alt_dic[chrom][pos] = [m1,m2,f1,f2,str(2)]
				elif f2 == m1:
					alt = f1
					if snp1 == alt or snp2 == alt:
						alt_dic[chrom][pos] = [m1,m2,f1,f2,str(1)]
				elif f1 != m1 and f2 != m1:
					alt = f1
					if snp1 == alt or snp2 == alt:
						alt_dic[chrom][pos] = [m1,m2,f1,f2,str(1)]
			elif m2 != ".":
				if f1 != m1 and f1 != m2:
					alt = f1
					if snp1 == alt or snp2 == alt:
						alt_dic[chrom][pos] = [m1,m2,f1,f2,str(1)]
				elif f2 != m1 and f2 != m2:
					alt = f2
					if snp1 == alt or snp2 == alt:
						alt_dic[chrom][pos] = [m1,m2,f1,f2,str(2)]
				elif m1 != f1 and m1 != f2 and m2 != f1 and m2 != f2:
					alt = f1
					if snp1 == alt or snp2 == alt:
						alt_dic[chrom][pos] = [m1,m2,f1,f2,str(1)]
			else:
				alt = "."

#order the dictionary 
oaltdic = alt_dic
for g in alt_dic:
	oaltdic[g] = OrderedDict(sorted(alt_dic[g].items()))
oaltdic = OrderedDict(sorted(oaltdic.items()))

#write dictionary in output file
for ch in oaltdic:
	for p in oaltdic[ch]:
		M1 = oaltdic[ch][p][0]
		M2 = oaltdic[ch][p][1]
		F1 = oaltdic[ch][p][2]
		F2 = oaltdic[ch][p][3]
		Sample = oaltdic[ch][p][4]
		text = ch+"\t"+str(p)+"\t"+M1+"\t"+M2+"\t"+F1+"\t"+F2+"\t"+Sample

		o.write(text+"\n")

afile.close()
o.close()
		















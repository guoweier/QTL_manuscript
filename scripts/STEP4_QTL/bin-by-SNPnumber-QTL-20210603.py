from __future__ import division
import os, sys, math, time
from optparse import OptionParser
from collections import OrderedDict

# This script bin SNPs for QTL analysis. 
# The SNP markers are binned by SNP NUMBER along the chromosomes for each sample. 
# It can be used for removing "noisy" SNPs, which are the sparse ones at opposite haplotype of the majority. 

# Logic: 
# a. By setting the number of SNPs in each bin, calculate the frequency of haplotypes among SNPs per bin. 
# b. If it satifies the setting threshold (90/95%), the higher freq haplotype is recorded as the bin-haplotype. 
# c. If it does not satify the setting threshold, it will show an "NA" in the output list. 

# Input: step3 ouput file from genotyping pipeline, with both reference and alternative alleles recorded. (script designed 20210602)

# Output: 
# Chrom(Chr01), Start(1), End(50), Haplotype-Sample1(1/2/NA), Freq-H1-Sample1(95%), Freq-H2-Sample1(5%), more samples

# Update 20210708: Record SNPs genomic positions, for calculating bin's middle locations (show in Mb) 
# This is used to make physical map in QTL analysis. 

usage = "\n\n%prog"
parser = OptionParser(usage=usage)

parser.add_option("-a", "--allele", dest="a", help="Input step3-both-allele-record file")
parser.add_option("-t", "--threshold", dest="t", type="float", default=0.9, help="Threshold for calling one haplotype")
parser.add_option("-s", "--size", dest="s", type="int", default=50, help="The number of SNPs in each bin (bin size)")
#parser.add_option("-f", "--female", dest="f", type="float", default=1.0, help="whether the input marker list is male or female haplotypes. Male=1, Female=2. Default is male")
parser.add_option("-o", "--outputfile", dest="o", help="Output file.")

(opt, args) = parser.parse_args()

#input parameters
afile = open(opt.a)
thresh = opt.t
binna = opt.s
o = open(opt.o,"w")

Libs = {}
Bins = {}
Pos = {}
Chrom = []
count = 0

#read input file header
#record sample name
#set dictionary for all samples
line = afile.readline()
line = line.split("\n")[0]
l = line.split("\t")
for i in l:
	if i.split("-")[0] == "Haplotype":
        	lib = i[10:]
        	Libs[lib,count] = {}
    	count += 1

#set bins
for line in afile:
	line = line.split("\n")[0]
	ln = line.split("\t")
	chrom = ln[0]
	pos = int(ln[1])
	#set SNP index
	if chrom not in Chrom:
		SNP = 1
		Chrom.append(chrom)
	SNP += 1
	bins = SNP//binna
	if len(str(bins)) < 7:
		binS = '0'*(7-len(str(bins)))+str(bins)
	else:
		binS = str(bins)
	BinNum = int(chrom.split("Chr")[1]+str(binS))
	if BinNum not in Bins:
		Bins[BinNum] = [chrom,bins*binna+1,(bins*binna)+binna]
		for Lib in Libs:
            Libs[Lib][BinNum] = {}
			Libs[Lib][BinNum]["H1"] = 0
            Libs[Lib][BinNum]["H2"] = 0
            Libs[Lib][BinNum]["ND"] = 0
    	# update 20210708: record SNPs genomic positions
    	if BinNum not in Pos:
    		Pos[BinNum] = []
    	Pos[BinNum] += [pos]
    	for i in Libs:
    		if ln[i[1]] == "1":
    			Libs[i][BinNum]["H1"] += 1
    		elif ln[i[1]] == "2":
    			Libs[i][BinNum]["H2"] += 1
    		else:
    			Libs[i][BinNum]["ND"] += 1

#calculate bin haplotypes
for lib in Libs:
	for binum in Libs[lib]:
		h1 = Libs[lib][binum]["H1"]
		h2 = Libs[lib][binum]["H2"]
		nd = Libs[lib][binum]["ND"]
		h1freq = h1/(h1+h2+nd)
		h2freq = h2/(h1+h2+nd)
		if (h1 > h2) and (h1freq >= thresh): 
			Libs[lib][binum]["Haplotype"] = "1"
			Libs[lib][binum]["Freq-H1"] = round(h1freq,3)
			Libs[lib][binum]["Freq-H2"] = round(h2freq,3)
			Libs[lib][binum]["EffectSNP"] = h1+h2
		elif (h1 < h2) and (h2freq >= thresh):
			Libs[lib][binum]["Haplotype"] = "2"
			Libs[lib][binum]["Freq-H1"] = round(h1freq,3)
			Libs[lib][binum]["Freq-H2"] = round(h2freq,3)
			Libs[lib][binum]["EffectSNP"] = h1+h2
		else:
			Libs[lib][binum]["Haplotype"] = "."
			Libs[lib][binum]["Freq-H1"] = round(h1freq,3)
			Libs[lib][binum]["Freq-H2"] = round(h2freq,3)
			Libs[lib][binum]["EffectSNP"] = h1+h2

#calculate bin's middle location (update 20210708)
for item in Pos:
	midpos = round(((min(Pos[item]) + max(Pos[item]))/2/1000000),3)
	Bins[item] += [midpos]

afile.close()

#sort dictionary
sortedbins = Bins.keys()
sortedbins.sort()
sortedlibs = Libs.keys()
sortedlibs.sort()
sortedpos = Pos.keys()
sortedpos.sort()

#set output file header
header = "Chrom\tStart\tEnd\tPos"
for Sample in sortedlibs:
    header += "\t"+"Haplotype-"+Sample[0]
    header += "\t"+"Freq-H1-"+Sample[0]
    header += "\t"+"Freq-H2-"+Sample[0]
    header += "\t"+"Effect-SNP-"+Sample[0]
o.write(header+'\n')

#write data into output file
for i in sortedbins:
    Chrom = str(Bins[i][0])
    Start = str(Bins[i][1])
    End = str(Bins[i][2])
    Pos = str(Bins[i][3])
    text = Chrom + '\t'+ Start + '\t' + End + '\t' + Pos
    
    for Sample in sortedlibs:
        haplo = Libs[Sample][i]["Haplotype"]
        freq1 = Libs[Sample][i]["Freq-H1"]
        freq2 = Libs[Sample][i]["Freq-H2"]
        effSNP = Libs[Sample][i]["EffectSNP"]
        text += "\t"+str(haplo)+"\t"+str(freq1)+"\t"+str(freq2)+"\t"+str(effSNP)

    o.write(text+'\n')

o.close()



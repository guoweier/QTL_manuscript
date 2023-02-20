import os, sys, math
from optparse import OptionParser
import threading, multiprocessing
from multiprocessing import Process, Queue

usage = "\n\n%prog"
parser = OptionParser(usage=usage)

parser.add_option("-a", "--allele", dest="a", help="Input step3-genotyping output file")
parser.add_option("-b", "--bin", dest="b", type="int", help="the bin size")
parser.add_option("-s", "--snp", dest="s", type="int", help="the least number of SNPs in each bin")
parser.add_option("-o", "--outputfile", dest="o", help="Output file.")

(opt, args) = parser.parse_args()

#open files
afile = open(opt.a)
binna = opt.b
snum = opt.s
o = open(opt.o, "w")

count = 0
validbin = 0
Libs = {}
Bins = {}
SNPs = {}

ahead = afile.readline()
ahead = ahead.split("\n")[0]
ahd = ahead.split("\t")
for item in ahd:
	if item.split("-")[0] == "TotalCov":
		lib = str(item[9:])
		Libs[lib,count] = {}
	count += 1
print(Libs)

for line in afile:
	line = line.split("\n")[0]
	l = line.split("\t")
	chrom = l[0]
	pos = int(l[1])
	bins = pos/binna

	if len(str(bins)) < 7:
		binS = '0'*(7-len(str(bins)))+str(bins)
	else:
		binS = str(bins)
	BinNum = int(chrom.split("Chr")[1]+str(binS))
	if BinNum not in Bins:
		Bins[BinNum] = [chrom,bins*binna,(bins*binna)+binna]
		for Lib in Libs:
            Libs[Lib][BinNum] = {}
            Libs[Lib][BinNum]['Cov'] = 0
            Libs[Lib][BinNum]['ExpBCov'] = 0
            Libs[Lib][BinNum]['ObsBCov'] = 0
        SNPs[BinNum] = 0
    for i in Libs:
        if l[i[1]+2] != '.':
            # add the cov
            Libs[i][BinNum]['Cov'] += int(l[i[1]])
            # add the expected number of instance of the male allele (genotype x cov)
            Libs[i][BinNum]['ExpBCov'] += int(float(l[i[1]])*(1-float(l[i[1]+1])))
            # add the observed number of instances of the B allele
            Libs[i][BinNum]['ObsBCov'] += int(float(l[i[1]])-float(l[i[1]+2]))
            SNPs[BinNum] += 1

for item in SNPs:
    if SNPs[item] < snum:
        for lib in Libs:
            Libs[lib] = {key:val for key, val in Libs[lib].items() if key != item}



afile.close()

sortedbins = Bins.keys()
sortedbins.sort()

sortedlibs = Libs.keys()
sortedlibs.sort()

header = "Chrom\tStart\tEnd"

for Sample in sortedlibs:
    header += "\t"+Sample[0]+"-Cov"

for Sample in sortedlibs:
    header += "\t"+Sample[0]+"-ObsMale"

for Sample in sortedlibs:
    header += "\t"+Sample[0]+"-CalcMale"

o.write(header+'\n')

for i in sortedbins:
    Chrom = str(Bins[i][0])
    Start = str(Bins[i][1])
    End = str(Bins[i][2])
    text = Chrom + '\t'+ Start + '\t' + End
    
    for Sample in sortedlibs:
        coverage = Libs[Sample][i]['Cov']
        text += '\t'+str(coverage)
    
    for Sample in sortedlibs:
        coverage = Libs[Sample][i]['Cov']
        if coverage > 0:
            ObsPerB = 100*(Libs[Sample][i]['ObsBCov'])/(Libs[Sample][i]['Cov'])
        else:
            ObsPerB =  "."
        text += '\t'+str(ObsPerB)

    for Sample in sortedlibs:
        coverage = Libs[Sample][i]['Cov']
        if coverage > 0:
            CalcPerB = 100*(Libs[Sample][i]['ExpBCov'])/(Libs[Sample][i]['Cov'])
        else:
            CalcPerB =  "."
        text += '\t'+str(CalcPerB)

    o.write(text+'\n')

o.close()






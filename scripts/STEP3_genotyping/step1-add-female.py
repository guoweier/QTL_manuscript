import os, sys, math, time
from optparse import OptionParser

# This script is going to add another parent marker list to one parent haplotype list based on positions.
# If the parent haplotype list is male, then add female for those SNP positions; if the parent haplotype list is female, then add male for those SNP positions.

# Input:
# 1. Parent haplotype list:
#   a. separated-haplotypes-list output (generated from step3 in build-haplotype) need to combine all blocks into a complete haplotype for every chromosome. This step needs to be done manually.
#   b. Format: Chrom+Pos+M1+M2 (Male); Chrom+Pos+F1+F2 (Female)
# 2. Original marker list before build-haplotype, which containing both male and female markers for all positions. (Add female: put-in selected-male; Add male: put-in selected-female)

# Output: Chrom+Pos+M1+M2+F1+F2

usage = "\n\n%prog"
parser = OptionParser(usage=usage)

parser.add_option("-s", "--snp", dest="s", help="Input marker list file")
parser.add_option("-f", "--female", dest="f", type="float", default=1.0, help="whether the input marker list is male or female haplotypes. Male=1, Female=2. Default is male")
parser.add_option("-m", "--marker", dest="m", help="Input selected marker file (containing both male and female markers)")
parser.add_option("-o", "--outputfile", dest="o", help="Output file.")

(opt, args) = parser.parse_args()

# open all files
marker = open(opt.m)
SNPfile = open(opt.s)
o = open(opt.o, "w")
mf = opt.f

SNPs = {}

SNPhead = SNPfile.readline()

for mine in SNPfile:
    mine = mine[:-1]
    m = mine.split("\t")
    if m[0] == m[1] == "-":
        continue
    else:
        SNPChrom = m[0]
        SNPPos = m[1]
        SNP1 = m[2]
        SNP2 = m[3][0]
        SNPMega = int(SNPPos)/1000000
        if SNPChrom not in SNPs:
            SNPs[SNPChrom] = {}
        if SNPMega not in SNPs[SNPChrom]:
            SNPs[SNPChrom][SNPMega] = {}
        SNPs[SNPChrom][SNPMega][SNPPos]=[SNP1,SNP2]

SNPfile.close()

# go check for selected female and male marker list
line = marker.readline()
l = line.split("\t")
if line[0] == "C":
    count = 0
    header = "Chrom\tPos\tM1\tM2\tF1\tF2"
    o.write(header+"\n")

for line in marker:
    line = line[:-1]
    l = line.split("\t")
    chrom = l[0]
    pos = l[1]
    mega = int(pos)/1000000
    text = ""
    if chrom in SNPs:
        if mega in SNPs[chrom]:
            if pos in SNPs[chrom][mega]:
                mf1 = SNPs[chrom][mega][pos][0]
                mf2 = SNPs[chrom][mega][pos][1]
                if mf = 1:
                    text = chrom+"\t"+str(pos)+"\t"+mf1+"\t"+mf2+"\t"+l[4]+"\t"+l[5]
                    o.write(text+"\n")
                elif mf = 2:
                    text = chrom+"\t"+str(pos)+"\t"+l[2]+"\t"+l[3]+"\t"+mf1+"\t"+mf2
                    o.write(text+"\n")

o.close()
marker.close()

import os, sys, math, time
from optparse import OptionParser

# This script uses parent marker haplotype list to genotyping every RNA sequencing sample.

# Input:
# 1. allele file: output file from step2-haplotype-call-allele.
#   a. Format: Chrom+Pos+M1+M2+F1+F2+Sample-info(SNPtype+SNP1+SNP2+TotalCov+%M1+CovM1)

# Output:
# 1. Chrom+Pos+M1+M2+F1+F2
# 2. Every sample genotype for that SNP marker. 

usage = "\n\n%prog"
parser = OptionParser(usage=usage)

parser.add_option("-a", "--allele", dest="a", help="Input step2-haplotype-call-allele output file")
parser.add_option("-f", "--female", dest="f", type="float", default=1.0, help="whether the input marker list is male or female haplotypes. Male=1, Female=2. Default is male")
parser.add_option("-o", "--outputfile", dest="o", help="Output file.")

(opt, args) = parser.parse_args()

# open files
allelefile = open(opt.a)
o = open(opt.o,"w")
mf = opt.f

indices = []
count = 0

line = allelefile.readline()
l = line.split("\t")
header = "Chrom\tPos\tM1\tM2\tF1\tF2"
for i in l:
    if i.split("-")[0] == "TotalCov":
        indices.append(count)
        lib = i[9:]
        header += "\t"+lib
    count += 1
o.write(header+"\n")

for line in allelefile:
    l = line.split("\t")
    chrom = l[0]
    pos = l[1]
    m1 = l[2]
    m2 = l[3]
    f1 = l[4]
    f2 = l[5]
    text = chrom+"\t"+pos+"\t"+m1+"\t"+m2+"\t"+f1+"\t"+f2
    # for male haplotype
    if mf == 1:
        for i in indices:
            snp1 = l[i-2]
            snp2 = l[i-1]
            if snp2 == ".":
                SnpM = snp1
            elif snp2 != ".":
                if f2 == ".":
                    if snp1 == f1:
                        SnpM = snp2
                    elif snp2 == f1:
                        SnpM = snp1
                    else:
                        SnpM = "."
                elif f2 != ".":
                    if ((snp1 == m1 or snp1 == m2) and (snp1 != f1 and snp1 != f2)) or ((snp2 == f1 or snp2 == f2) and (snp2 != m1 and snp2 != m2)):
                        SnpM = snp1
                    elif ((snp1 == f1 or snp1 == f2) and (snp1 != m1 and snp2 != m2)) or ((snp2 == m1 or snp2 == m2) and (snp2 != f1 and snp2 != f2)):
                        SnpM = snp2
                    else:
                        SnpM = "."
            text += "\t"+SnpM
        o.write(text+"\n")
    # for female haplotype
    elif mf == 2:
        for i in indices:
            snp1 = l[i-2]
            snp2 = l[i-1]
            if snp2 == ".":
                SnpF = snp1
            elif snp2 != ".":
                if m2 == ".":
                    if snp1 == m1:
                        SnpF = snp2
                    elif snp2 == m1:
                        SnpF = snp1
                    else:
                        SnpF = "."
                elif m2 != ".":
                    if ((snp1 == f1 or snp1 == f2) and (snp1 != m1 and snp1 != m2)) or ((snp2 == m1 or snp2 == m2) and (snp2 != f1 and snp2 != f2)):
                        SnpF = snp1
                    elif ((snp1 == m1 or snp1 == m2) and (snp1 != f1 and snp2 != f2)) or ((snp2 == f1 or snp2 == f2) and (snp2 != m1 and snp2 != m2)):
                        SnpF = snp2
                    else:
                        SnpF = "."
            text += "\t"+SnpF
        o.write(text+"\n")

# close files
o.close()
allelefile.close()

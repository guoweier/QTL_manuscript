import os, sys, math, time
from optparse import OptionParser


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
line = line.split("\n")[0]
l = line.split("\t")
header = "Chrom\tPos\tM1\tM2\tF1\tF2"
for i in l:
    if i.split("-")[0] == "TotalCov":
        indices.append(count)
        lib = i[9:]
        if mf == 1:
            libsnp = "Male-"+lib
            tcov = "TotalCov-"+lib
            hapsum = "Haplotype-"+lib
        elif mf == 2:
            libsnp = "Female-"+lib
            tcov = "TotalCov-"+lib
            hapsum = "Haplotype-"+lib
        header += "\t"+libsnp+"\t"+tcov+"\t"+hapsum
    count += 1
o.write(header+"\n")

for line in allelefile:
    line = line.split("\n")[0]
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
            totalcov = l[i]
            if snp1 == ".":
                SnpM = "."
                hap = "."
            elif snp2 == ".":
                if snp1 == f1:
                    if snp1 == m1:
                        SnpM = snp1
                        hap = str(1)
                    elif snp1 == m2:
                        SnpM = snp1
                        hap = str(2)
                    else:
                    	SnpM = "."
                        hap = "."
                elif snp1 != f1:
                    if snp1 == m1:
                        SnpM = snp1
                        hap = str(1)
                    elif snp1 == m2:
                        SnpM = snp1
                        hap = str(2)
                    else:
                        SnpM = "."
                        hap = "."
            elif snp2 != ".":
                if f2 == ".":
                    if snp1 == f1:
                        SnpM = snp2
                        if snp2 == m1:
                            hap = str(1)
                        elif snp2 == m2:
                            hap = str(2)
                    elif snp2 == f1:
                        SnpM = snp1
                        if snp1 == m1:
                            hap = str(1)
                        elif snp1 == m2:
                            hap = str(2)
                    else:
                        SnpM = "."
                        hap = "."
                elif f2 != ".":
                    if ((snp1 == m1 or snp1 == m2) and (snp1 != f1 and snp1 != f2)) or ((snp2 == f1 or snp2 == f2) and (snp2 != m1 and snp2 != m2)):
                        SnpM = snp1
                        if snp1 == m1:
                            hap = str(1)
                        elif snp1 == m2:
                            hap = str(2)
                    elif ((snp1 == f1 or snp1 == f2) and (snp1 != m1 and snp2 != m2)) or ((snp2 == m1 or snp2 == m2) and (snp2 != f1 and snp2 != f2)):
                        SnpM = snp2
                        if snp2 == m1:
                            hap = str(1)
                        elif snp2 == m2:
                            hap = str(2)
                    else:
                        SnpM = "."
                        hap = "."
            text += "\t"+SnpM+"\t"+totalcov+"\t"+hap
        o.write(text+"\n")
    # for female haplotype
    elif mf == 2:
        for i in indices:
            snp1 = l[i-2]
            snp2 = l[i-1]
            totalcov = l[i]
            if snp1 == ".":
                SnpF = "."
                hap = "."
            elif snp2 == ".":
                if snp1 == m1:
                    if snp1 == f1:
                        SnpF = snp1
                        hap = str(1)
                    elif snp1 == f2:
                        SnpF = snp1
                        hap = str(2)
                    else:
                    	SnpF = "."
                        hap = "."
                elif snp1 != m1:
                    if snp1 == f1:
                        SnpF = snp1
                        hap = str(1)
                    elif snp1 == f2:
                        SnpF = snp1
                        hap = str(2)
                    else:
                        SnpF = "."
                        hap = "."
            elif snp2 != ".":
                if m2 == ".":
                    if snp1 == m1:
                        SnpF = snp2
                        if snp2 == f1:
                            hap = str(1)
                        elif snp2 == f2:
                            hap = str(2)
                    elif snp2 == m1:
                        SnpF = snp1
                        if snp1 == f1:
                            hap = str(1)
                        elif snp1 == f2:
                            hap = str(2)
                    else:
                        SnpF = "."
                        hap = "."
                elif m2 != ".":
                    if ((snp1 == f1 or snp1 == f2) and (snp1 != m1 and snp1 != m2)) or ((snp2 == m1 or snp2 == m2) and (snp2 != f1 and snp2 != f2)):
                        SnpF = snp1
                        hap = str(1)
                    elif ((snp1 == m1 or snp1 == m2) and (snp1 != f1 and snp2 != f2)) or ((snp2 == f1 or snp2 == f2) and (snp2 != m1 and snp2 != m2)):
                        SnpF = snp2
                        hap = str(2)
                    else:
                        SnpF = "."
                        hap = "."
            text += "\t"+SnpF+"\t"+totalcov+"\t"+hap
        o.write(text+"\n")

# close files
o.close()
allelefile.close()
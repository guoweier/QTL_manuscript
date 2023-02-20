import os, sys, math, time
from optparse import OptionParser

usage = "\n\n%prog"
parser = OptionParser(usage=usage)

parser.add_option("-H", "--Haplotype", dest="haplo", help="Input haplotype file.")
#parser.add_option("-f", "--outputfolder", dest="f", help="Output main folder name.")
parser.add_option("-o", "--outputfile", dest="o", help="Output file.")


(opt, args) = parser.parse_args()

#files
haplofile = open(opt.haplo)
o = open(opt.o,'w')
chroms = []
haplotb = []
outputtb = []

haplohead = haplofile.readline()

#adjust input file
for line in haplofile:
    ln = line.split('\t')
    chrom = ln[0]
    pos = ln[1]
    m1 = ln[2]
    m2 = ln[3]
    if chrom not in chroms:
        chroms.append(chrom)
    item = [chrom, pos, m1, m2]
    haplotb.append(item)

haplofile.close()

#write output file head
header = haplohead
o.write(header)

for i in range(len(haplotb)):
    if haplotb[i][0] == haplotb[i][1] == "-":
        txt = haplotb[i][0]+'\t'+haplotb[i][1]+'\t'+haplotb[i][2]+'\t'+haplotb[i][3][0]
        outputtb.append(txt)
    else:
        Ch1 = haplotb[i][0]
        Pos1 = int(haplotb[i][1])
        Ma1 = haplotb[i][2]
        Ma2 = haplotb[i][3]
        ni = i+1
        if ni < len(haplotb) and haplotb[ni][0] != "-":
            Ch2 = haplotb[ni][0]
            Pos2 = int(haplotb[ni][1])
            Mb1 = haplotb[ni][2]
            Mb2 = haplotb[ni][3]
            if Ch1 == Ch2:
                text1 = Ch1+'\t'+str(Pos1)+'\t'+Ma1+'\t'+Ma2[0]
                text2 = Ch2+'\t'+str(Pos2)+'\t'+Mb1+'\t'+Mb2[0]
                if Pos2 > Pos1:
                    if text1 not in outputtb:
                        outputtb.append(text1)
                    outputtb.append(text2)
                elif Pos2 <= Pos1:
                    outputtb.append('----')
                    outputtb.append(text2)
            elif Ch1 != Ch2:
                outputtb.append('----\n----')

for m in outputtb:
    o.write(m+'\n')

o.close()

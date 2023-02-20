from __future__ import division
import os, sys, math, time
from optparse import OptionParser
from collections import OrderedDict

usage = "\n\n%prog"
parser = OptionParser(usage=usage)

parser.add_option("-s", "--selected", dest="s", help="Input selected male het haplotype file")
parser.add_option("-o", "--outputfile", dest="o", help="Output percentage file.")

(opt, args) = parser.parse_args()

#selected haplotype file
selectfile = open(opt.s)
#output file
o = open(opt.o, 'w')

indices = []
count = 0
haplo = [] #result haplotype collecting table
check = [] #haplotype collecting check table
selectdic = {}

#read header and record sample numbers
selecthead = selectfile.readline()
shead = selecthead.split('\t')
for item in shead:
    if item.split('-')[0] == 'Male':
        lib = item[5:]
        indices.append([lib, count])
    count += 1

#write header for output file
header = "Chrom\tPos\tM1\tM2\tPercentage"
o.write(header+'\n')

for line in selectfile:
    sample = []
    #read lines and record into dictory
    ln = line.split('\t')
    chrom = int(ln[0][3:])
    pos = int(ln[1])
    m1 = ln[2]
    m2 = ln[3]
    for a in indices:
        sample.append(ln[a[1]])
    if chrom not in selectdic:
        selectdic[chrom] = {}
    selectdic[chrom][pos] = [m1, m2, sample]
orderdic = selectdic
for g in selectdic:
    orderdic[g] = OrderedDict(sorted(selectdic[g].items()))
selectfile.close()

for c in orderdic:
    ch = "Chr"+str(c)
    print(c)
    markers = [] #set a table for markers position and het alleles
    for item in orderdic[c]:
        m1 = orderdic[c][item][0]
        m2 = orderdic[c][item][1]
        markers.append([item,m1,m2])
    for i in range(len(markers)):
        samplesnp = [] #first marker sample snps
        ni = i #prepare for next marker

        Pos = markers[i][0]
        m1 = markers[i][1]
        m2 = markers[i][2]
        for s in range(len(orderdic[c][Pos][2])):
            samplesnp.append(orderdic[c][Pos][2][s])

        nsamplesnp = [] #second marker sample snps
        ni = ni+1
        if ni < len(markers):
            totalsample = [] #sample snps combination collecting table
            nPos = markers[ni][0]
            nm1 = markers[ni][1]
            nm2 = markers[ni][2]
            com1 = m1+nm1
            com2 = m2+nm2
            recom1 = m1+nm2
            recom2 = m2+nm1
            norecom = 0
            recom = 0
            for ns in range(len(orderdic[c][nPos][2])):
                nsamplesnp.append(orderdic[c][nPos][2][ns])
            for a in range(len(samplesnp)):
                samplecomb = str(samplesnp[a])+str(nsamplesnp[a])
                totalsample.append(samplecomb)
            for snp in totalsample:
                if snp == com1 or snp == com2:
                    norecom += 1
                elif snp == recom1 or snp == recom2:
                    recom += 1
                else:
                    continue
            prob_norecom = float(norecom/len(totalsample))*100
            prob_recom = float(recom/len(totalsample))*100
            if prob_norecom >= prob_recom:
                snp1ch = [ch, str(Pos), m1, m2]
                snp1pch = [ch, str(Pos), m2, m1]
                snp2ch = [ch, str(nPos), nm1, nm2]
                snp2pch = [ch, str(nPos), nm2, nm1]
                snp1 = [ch, str(Pos), m1, m2, str(prob_norecom)]
                snp1p = [ch, str(Pos), m2, m1, str(prob_norecom)]
                snp2 = [ch, str(nPos), nm1, nm2, "."]
                snp2p = [ch, str(nPos), nm2, nm1, "."]
                if (snp1ch not in check) and (snp1pch not in check):
                    check.append(snp1ch)
                    check.append(snp2ch)
                    haplo.append(snp1)
                    haplo.append(snp2)
                elif snp1ch in check:
                    check.append(snp1ch)
                    check.append(snp2ch)
                    haplo.append(snp1)
                    haplo.append(snp2)
                elif snp1pch in check:
                    check.append(snp1pch)
                    check.append(snp2pch)
                    haplo.append(snp1p)
                    haplo.append(snp2p)
            elif prob_norecom < prob_recom:
                snp1ch = [ch, str(Pos), m1, m2]
                snp1pch = [ch, str(Pos), m2, m1]
                snp2ch = [ch, str(nPos), nm1, nm2]
                snp2pch = [ch, str(nPos), nm2, nm1]
                snp1 = [ch, str(Pos), m1, m2, str(prob_recom)]
                snp1p = [ch, str(Pos), m2, m1, str(prob_recom)]
                snp2 = [ch, str(nPos), nm1, nm2, "."]
                snp2p = [ch, str(nPos), nm2, nm1, "."]
                if (snp1ch not in check) and (snp1pch not in check):
                    check.append(snp1ch)
                    check.append(snp2pch)
                    haplo.append(snp1)
                    haplo.append(snp2p)
                elif snp1ch in check:
                    check.append(snp1ch)
                    check.append(snp2pch)
                    haplo.append(snp1)
                    haplo.append(snp2p)
                elif snp1pch in check:
                    check.append(snp1pch)
                    check.append(snp2ch)
                    haplo.append(snp1p)
                    haplo.append(snp2)

#write haplo table info into txt file
for unit in haplo:
    text = unit[0]+'\t'+unit[1]+'\t'+unit[2]+'\t'+unit[3]+'\t'+unit[4]
    o.write(text+'\n')

o.close()
selectfile.close()

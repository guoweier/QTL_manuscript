from __future__ import division
import os, sys, math, time
from optparse import OptionParser
from collections import OrderedDict

# This is the STEP2 of generating haplotypes based on snp markers (male).
# The logic of this script:
# 1. Combine every adjacent two markers. Marker1: A T; Marker2: A G.
# 2. Haplotype can be either AA/TG or AG/TA
# 3. Check alleles combination for these two markers on each sample.
# 4. Count samples allele combination for AA/TG and AG/TA
# 5. Calculate the percentage of two types. Either one of both is larger than set threshold (default 90%), that combine is the haplotype.
# 6. If both of them are smaller than threshold, skip the lower line. Check for the first line with the next following line.
# 7. If skipped "lower" lines beyond set times (default 5), skip the first line, and start with the first skipped line again.
# 8. All the skipped lines would be recorded in odd file.

# Input file:
#  Selected snp list. The list only contains heterozygous snps for male.
#  Format:
#   1. Chrom: Chromosome ref
#   2. Pos: SNP marker position
#   3. M1: The most frequent allele for that marker in P.nigra
#   4. M2: The second most frequent allele for that marker in P.nigra
#   5. F1: The most frequent allele for that marker in P.deltoide (not used here)
#   6. F2: The second most frequent allele for that marker in P.deltoide (not used here)
#   7. Samples info, include: allele from Male, allele from Female, Total reads coverage for that position.

# Output files:
#  1. The list of haplotypes combination.
#     Format: Chrom, Pos, M1, M2. (M1 and M2 columns are two haplotypes)
#  2. The list of odd lines that been skipped.
#     Format: Chrom, Pos, M1, M2, Inver(1=skipped as the first line; 2=skipped as the 2-5 lines)

usage = "\n\n%prog"
parser = OptionParser(usage=usage)

parser.add_option("-s", "--selected", dest="s", help="Input selected male het haplotype file")
parser.add_option("-t", "--threshold", dest="t", type = "float", default = 90.0, help="minimun percentage for determine linkage")
parser.add_option("-n", "--nextline", dest="next", type = "int", default = 5, help="Maximum number of next line skipped")
parser.add_option("-o", "--outputfile", dest="o", help="Output file.")
parser.add_option("-d", "--oddname", dest="odd", help = "File name for odd lines.")

(opt, args) = parser.parse_args()

#selected haplotype file
selectfile = open(opt.s)
#output file
o = open(opt.o, 'w')
oddf = open(opt.odd, 'w')
#parameters put-in
thresh = opt.t
nextnum = opt.next

indices = []
count = 0
haplo = [] #result haplotype collecting table
odd = [] #table collecting weird markers that the combination less than thresh (90%)
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
header = "Chrom\tPos\tM1\tM2"
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
        j = 0 #parameter for checking insersion
        ni = i #prepare for next marker
        fin = False #judge whether need to skip the line

        Pos = markers[i][0]
        m1 = markers[i][1]
        m2 = markers[i][2]
        for s in range(len(orderdic[c][Pos][2])):
            samplesnp.append(orderdic[c][Pos][2][s])

        while bool(fin) == False and j < nextnum:
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
                if float(norecom/len(totalsample))*100 >= float(thresh):
                    snp1 = [ch, str(Pos), m1, m2]
                    snp1p = [ch, str(Pos), m2, m1]
                    snp2 = [ch, str(nPos), nm1, nm2]
                    snp2p = [ch, str(nPos), nm2, nm1]
                    if (snp1 not in haplo) and (snp1p not in haplo):
                        haplo.append(snp1)
		    elif ((snp1 != haplo[(len(haplo)-1)]) and (snp1p != haplo[(len(haplo)-1)])):
                        if snp1 in haplo:
                            haplo.append(snp1)
                        elif snp1p in haplo:
                            haplo.append(snp1p)
                    if snp1 in haplo:
                        haplo.append(snp2)
                    elif snp1p in haplo:
                        haplo.append(snp2p)
                    fin = True
                elif float(recom/len(totalsample))*100 >= float(thresh):
                    snp1 = [ch, str(Pos), m1, m2]
                    snp1p = [ch, str(Pos), m2, m1]
                    snp2 = [ch, str(nPos), nm1, nm2]
                    snp2p = [ch, str(nPos), nm2, nm1]
                    if (snp1 not in haplo) and (snp1p not in haplo):
                        haplo.append(snp1)
                    elif ((snp1 != haplo[(len(haplo)-1)]) and (snp1p != haplo[(len(haplo)-1)])):
                        if snp1 in haplo:
                            haplo.append(snp1)
                        elif snp1p in haplo:
                            haplo.append(snp1p)
                    if snp1 in haplo:
                        haplo.append(snp2p)
                    elif snp1p in haplo:
                        haplo.append(snp2)
                    fin = True
                elif float(norecom/len(totalsample))*100 < float(thresh) and float(recom/len(totalsample))*100 < float(thresh):
                    snp2 = [ch, str(nPos), nm1, nm2, str(2)]
                    odd.append(snp2)
                    fin = False
            j = j+1
        if fin == False:
            snp1 = [ch, str(Pos), m1, m2, str(1)]
            odd.append(snp1)
	    if len(haplo) == 0:
		continue
	    elif haplo[len(haplo)-1] == [ch, str(Pos), m1, m2] or haplo[len(haplo)-1] == [ch, str(Pos), m2, m1]:
		haplo.append(["-","-","-","-"])
        elif fin == True:
            continue

#write haplo table info into txt file
for unit in haplo:
    text = unit[0]+'\t'+unit[1]+'\t'+unit[2]+'\t'+unit[3]
    o.write(text+'\n')

#write odd table info into txt file
oddhead = "Chrom\tPos\tM1\tM2\tInver"
oddf.write(oddhead+'\n')
for unit in odd:
    textodd = unit[0]+'\t'+unit[1]+'\t'+unit[2]+'\t'+unit[3]+'\t'+unit[4]
    oddf.write(textodd+'\n')

#close files
o.close()
oddf.close()

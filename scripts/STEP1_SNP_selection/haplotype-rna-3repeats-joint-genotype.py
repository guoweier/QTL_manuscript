#! /usr/bin/env python

import os, sys, math
from optparse import OptionParser

usage=""
parser = OptionParser(usage=usage)
parser.add_option('-f', '--file', dest='f', help='Input txt for containing genotype information.')
parser.add_option('-r', '--rfile', dest='r', help='Input txt again.')

(opt, args) = parser.parse_args()

indexf = opt.f
rf = opt.r
f = open(indexf,'r')
r = open(rf, 'r')

# extract all the mRNA files into a list #
li1 = os.listdir('/share/comailab/sequencing/190118_J00113_0437_AH2KM3BBXY_run612_2019-01-18_H1667_1678P_Tsai/prep/output/')
li2 = os.listdir('/cato2pool/pop-mRNA-20lanes/20180808/combined/fqs/')
li3 = os.listdir('/cato2pool/pop-mRNA-20lanes/20181024/fqs/')
li4 = os.listdir('/cato2pool/readsets/NovaSeq_2019-03-22_POP_mRNASeq_20190225_HTsai/prepped/trimmed/')

path1 = '/share/comailab/sequencing/190118_J00113_0437_AH2KM3BBXY_run612_2019-01-18_H1667_1678P_Tsai/prep/output'
path2 = '/cato2pool/pop-mRNA-20lanes/20180808/combined/fqs'
path3 = '/cato2pool/pop-mRNA-20lanes/20181024/fqs'
path4 = '/cato2pool/readsets/NovaSeq_2019-03-22_POP_mRNASeq_20190225_HTsai/prepped/trimmed'
lsfilt={path1:[],path2:[],path3:[],path4:[]}

for fname in li1:
    if 'POP_mRNA' in fname:
        lsfilt[path1] += [fname]
for fname in li2:
    if 'POP_mRNA' in fname:
        lsfilt[path2] += [fname]
for fname in li3:
    if 'POP_mRNA' in fname:
        lsfilt[path3] += [fname]
for fname in li4:
    if 'POP_mRNA' in fname:
        lsfilt[path4] += [fname]
                

# set a unique list of genotypes #
genoli = []
for ln in f:
    ls = ln.split('\t')
    if ls[0] == "RNA Seq Set":
        continue
    elif ls[11] not in genoli:
        genoli.append(ls[11])
    elif ls[11] in genoli:
        continue


# write different fq files into one same fq file #
x=[]
for ln in r:
    if ln == "":
        break
    ls = ln.split('\t')
    if ls[0] == "RNA Seq Set":
        continue
    x.append([ls[1],ls[11]])
for item in genoli:
    print(item)
    fqli=[]
    fqfinal=open(item + '.fq', 'w')
    for a in x:
        if a[1] == item:
            print(a[0])
            fqli.append(a[0])
    for b in fqli:
        for paths in lsfilt:
            for fname in lsfilt[paths]:
                if b in fname:
                    os.chdir(paths)
                    fq1=open(fname)
                    while 1:
                        name = fq1.readline()
                        if name == "":
                            break
                        seq = fq1.readline()
                        plus = fq1.readline()
                        qual = fq1.readline()

                        fqfinal.write(name+seq+plus+qual)
    fqfinal.close()
    os.chdir('/cato2pool/weier-poplar/haplotype-project/rnaseq')                           



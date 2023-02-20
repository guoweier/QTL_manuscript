import os, sys, math, time
from optparse import OptionParser

# This is script is similar as call-alleles from comailab script-share.
# It uses SNP markers haplotype to call alleles in each samples from parsed-mpileup files.
# Input:
# 1. SNP haplotype list: Must be M1+M2+F1+F2. For male, M1 and M2 must be the haplotype ordering list. (See build-haplotype folder for details about how to get haplotype for parents) F1 and F2 are just female added into male haplotype list. For female, F1 and F2 must be the haplotype ordering list. M1 and M2 are just male added into female haplotype list. Use step1-add-other-parent in advance, then the SNP input format would be: Chrom+Pos+M1+M1+F1+F2
# 2. parsed-mpileup file: same parsed-mpileup file used in previous call-alleles in parents genotyping.

# Output: Similar to the output of call-alleles from comailab script-share.
# 1. Chromosome, Allele position, allele1 in Male, allele2 in Male, allele1 from Female, allele2 from Female (being '.' when female is homo).
# 2. 6 columns per sample:
#    a. SNP type: 1=homo, 2=het
#    b. SNP1: most common allele
#    c. SNP2: second allele
#    d. Total Cov: total coverage at that position
#    e. %M1: 1=homo for allele in Male, 0.5=het, 0=homo for allele in Female
#    f. CovM1: coverage of the allele from Male

usage = "\n\n%prog"
parser = OptionParser(usage=usage)

parser.add_option("-s", "--snp", dest="s", help="Input haplotype SNP marker file")
parser.add_option("-m", "--mpileup", dest="m", help="Input mpileup file")
parser.add_option("-f", "--female", dest="f", type="float", default=1.0, help="whether the input marker list is male or female haplotypes. Male=1, Female=2. Default is male")
parser.add_option("-H", "--hetmin", dest="het", type = "float", default = 10.0, help="Maximum percentage of second allele to be considered het defualt is 10% (10.0)")
parser.add_option("-o", "--outputfile", dest="o", help="Output file.")

(opt, args) = parser.parse_args()

mpup = open(opt.m) #mpileup file
o = open(opt.o,'w') #output file
mf = opt.f

## make a list of SNP position and their base call in the two parental genotypes ##

SNPfile = open(opt.s) #snp file
SNPs = {}

SNPhead = SNPfile.readline()

for mine in SNPfile:
    m = mine.split("\t")
    if m[0] == m[1] == "-":
        continue
    else:
        SNPChrom = m[0]
        SNPPos = m[1]
        SNPm1 = m[2]
        SNPm2 = m[3]
        SNPf1 = m[4]
        SNPf2 = m[5][0]
        SNPMega = int(SNPPos)/1000000
        if SNPChrom not in SNPs:
            SNPs[SNPChrom] = {}
        if SNPMega not in SNPs[SNPChrom]:
            SNPs[SNPChrom][SNPMega] = {}
        SNPs[SNPChrom][SNPMega][SNPPos]=[SNPm1,SNPm2,SNPf1,SNPf2]

SNPfile.close()

## Going through the mpileup now ##

indices = []
totalcount = 0
odd1 = 0

line = mpup.readline()
l = line.split()
if line[0] == "C":
    count = 0
    header = "Chrom\tPos\tM1\tM2\tF1\tF2"
    for m in l:
        if m.split('-')[0] == 'Cov':
            indices.append(count)
            lib = m[4:]
            header += '\tSnptype-'+lib+'\tSNP1-'+lib+'\tSNP2-'+lib+'\tTotalCov-'+lib+'\t%M1-'+lib+'\tCovM1-'+lib
        count += 1
    o.write(header+'\n')

for line in mpup:
    l = line.split()
    if l[0][0] == "C":
        chrom = l[0]
        pos = l[1]
        mega = int(pos)/1000000
        text=""
        if chrom in SNPs:
            if mega in SNPs[chrom]:
                if pos in SNPs[chrom][mega]:
                    m1 = SNPs[chrom][mega][pos][0]
                    m2 = SNPs[chrom][mega][pos][1]
                    f1 = SNPs[chrom][mega][pos][2]
                    f2 = SNPs[chrom][mega][pos][3]
                    text = chrom+"\t"+str(pos)+"\t"+m1+"\t"+m2+"\t"+f1+"\t"+f2
                    print(text)

                    for i in indices:
                        Cov = l[i]
                        if Cov == '.':
                            text += 6 * '\t.'
                        else:
                            totalcount += 1
                            snp1 = l[i-3]
                            snp2 = l[i-2]
                            snp3 = l[i-1]
                            if snp3 != '.':
                                Snptype = 3
                            elif snp2 != '.':
                                Snptype = 2
                            else:
                                Snptype = 1
                            Snp1 = snp1.split("_")[0]
                            Cov1 = int(round(float(snp1.split("_")[1])*float(Cov)/100))
                            if snp2 == '.':
                                Snp2 = '.'
                            elif float(snp2.split('_')[1]) <= opt.het:
                                Snp2 = '.'
                            else:
                                Snp2 = snp2.split('_')[0]
                                Cov2 = int(round(float(snp2.split('_')[1])*float(Cov)/100))
                            # for male haplotype
                            if mf == 1:
                                if Snp2 == '.':
                                    if Snp1 == m1:
                                        All = 1
                                        CovM1 = Cov
                                    elif Snp1 == m2:
                                        All = 0
                                        CovM1 = 0
                                elif Snp2 != '.':
                                    if Snp1 == m1:
                                        All = 0.5
                                        CovM1 = Cov1
                                    elif Snp1 == m2:
                                        All = 0.5
                                        CovM1 = Cov2
                                else:
                                    All = '.'
                                    odd1 += 1
                                    CovM1 = '.'
                            # for female haplotype
                            elif mf == 2:
                                if Snp2 == ".":
                                    if Snp1 == f1:
                                        All = 1
                                        CovM1 = Cov
                                    elif Snp1 == f2:
                                        All = 0
                                        CovM1 = 0
                                elif Snp2 != ".":
                                    if Snp1 == f1:
                                        All = 0.5
                                        CovM1 = Cov1
                                    elif Snp1 == f2:
                                        All = 0.5
                                        CovM1 = Cov2
                                else:
                                    All = "."
                                    odd1 += 1
                                    CovM1 = "."
                            text += '\t'+str(Snptype)+'\t'+str(Snp1)+'\t'+str(Snp2)+'\t'+str(Cov)+'\t'+str(All)+'\t'+str(CovM1)
                    text += '\n'
                    o.write(text)

o.close()
mpup.close()

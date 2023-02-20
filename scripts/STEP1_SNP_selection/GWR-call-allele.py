import os, sys, math, time
from optparse import OptionParser

## This script takes parsed-mpileup file as the input. ##
## Using the heterozygous SNPs list generated from whole-genome-sequencing poplar parents (both het) as the SNP reference. ##
## Het SNPs list format: Chr, Pos, Male-allele (m_A or m_A/T), Male-Cov, Female-allele (f_A or f_A/T), Female-Cov ##
## The alleles between male and female can be in following format: ##
## 1. M_A F_A/T 
## 2. M_A F_T 
## 3. M_A F_T/A 
## 4. M_A/T F_A
## 5. M_T/A F_A
## 6. M_A/T F_A/C

## The logic of this script: ##
## 1. Only one allele is found, and that allele belongs to one of the parent, that allele is called. 
## 2. If two alleles are found, and the second allele percentage is < 10% (can be changed), only the most common allele is called.
## 3. If two alleles are found, and the second allele percentage is > 10% (can be changed), two alleles are called.
## 4. Other genotypes are assigned as 'NA'.

## Output:
## 1. Chromosome, Allele position, allele1 in Male, allele2 in Male (being '.' when male is homo), allele1 from Female, allele2 from Female (being '.' when female is homo).
## 2. 6 columns per sample:
##    a. SNP type: 1=homo, 2=het
##    b. SNP1: most common allele
##    c. SNP2: second allele
##    d. Total Cov: total coverage at that position
##    e. %M: 1=homo for allele in Male, 0.5=het, 0=homo for allele in Female
##    f. CovM: coverage of the allele from Male


usage = "\n\n%prog"
parser = OptionParser(usage=usage)

parser.add_option("-s", "--snp", dest="s", help="Input GWR-SNP File")
parser.add_option("-m", "--mpileup", dest="m", help="Input mpileup file")
parser.add_option("-H", "--hetmin", dest="het", type = "float", default = 10.0, help="Maximum percentage of second allele to be considered het defualt is 10% (10.0)")
parser.add_option("-o", "--outputfile", dest="o", help="Output file.")


(opt, args) = parser.parse_args()

mpup = open(opt.m) #mpileup file
o = open(opt.o,'w') #output file

## make a list of SNP position and their base call in the two parental genotypes ##

SNPfile = open(opt.s) #snp file
SNPs = {}

SNPhead = SNPfile.readline()

for mine in SNPfile:
    m = mine.split('\t')
    if '.-' in m[2] or '.-' in m[4]:
        continue
    else:
        SNPChrom = m[0]
        SNPPos = m[1]
        if len(m[2]) == 3:
            SNP_M1 = m[2].split('_')[1]
            SNP_M2 = '.'
        elif len(m[2]) > 3:
            SNP_M1 = m[2].split('_')[1].split('/')[0]
            SNP_M2 = m[2].split('_')[1].split('/')[1]
            if SNP_M2[:1] == '.+' and SNP_M2[3] != SNP_M1:
                SNP_M2 = SNP_M2[3]
        if len(m[4]) == 3:
            SNP_F1 = m[4].split('_')[1]
            SNP_F2 = '.'
        elif len(m[4]) > 3:
            SNP_F1 = m[4].split('_')[1].split('/')[0]
            SNP_F2 = m[4].split('_')[1].split('/')[1]
            if SNP_F2[:1] == '.+' and SNP_F2[3] != SNP_F1:
                SNP_F2 = SNP_F2[3]
        SNPMega = int(SNPPos)/1000000
        if SNPChrom not in SNPs:
            SNPs[SNPChrom] = {}
        if SNPMega not in SNPs[SNPChrom]:
            SNPs[SNPChrom][SNPMega] = {}
        SNPs[SNPChrom][SNPMega][SNPPos]=[SNP_M1,SNP_M2,SNP_F1,SNP_F2]
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
            header += '\tSnptype-'+lib+'\tSNP1-'+lib+'\tSNP2-'+lib+'\tTotalCov-'+lib+'\t%M-'+lib+'\tCovM-'+lib
        count += 1
    o.write(header+'\n')

for line in mpup:
    l = line.split()
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
                text = chrom+'\t'+str(pos)+'\t'+m1+'\t'+m2+'\t'+f1+'\t'+f2

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
                        Snp1 = snp1.split('_')[0]
                        Cov1 = int(round(float(snp1.split("_")[1])*float(Cov)/100))
                        if snp2 == '.':
                            Snp2 = '.'
                        elif float(snp2.split('_')[1]) <= opt.het:
                            Snp2 = '.'
                        else:
                            Snp2 = snp2.split('_')[0]
                            Cov2 = int(round(float(snp2.split('_')[1])*float(Cov)/100))
                        if Snp2 == '.':
                            if Snp1 == m1 or Snp1 == m2:
                                All = 1
                                CovM = Cov
                            elif Snp1 == f1 or Snp1 == f2:
                                All = 0
                                CovM = 0
                        elif Snp2 != '.':
                            if (Snp1 == m1 or Snp1 == m2) and (Snp2 == f1 or Snp2 == f2):
                                All = 0.5
                                CovM = Cov1
                            elif (Snp1 == f1 or Snp1 == f2) and (Snp2 ==m1 or Snp2 == m2):
                                All =0.5
                                CovM =Cov2
                        else:
                            All = '.'
                            odd1 += 1
                            CovM = '.'
                        text += '\t'+str(Snptype)+'\t'+str(Snp1)+'\t'+str(Snp2)+'\t'+str(Cov)+'\t'+str(All)+'\t'+str(CovM)
                text += '\n'
                o.write(text)

o.close()
mpup.close()

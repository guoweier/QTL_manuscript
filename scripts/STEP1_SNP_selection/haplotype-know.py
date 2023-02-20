import os, sys, math, time
from optparse import OptionParser

## This script takes call-allele file as the input.
## Using whole-genome-sequencing SNPs as the reference, set the haplotypes for each genotyping sample.

## The logic of this script:
## 1. If the sample is homozygous at SNP position, it inherited same allele from paternal and maternal.
## 2. If the sample is heterozygous at SNP position, in inherited different allele from paternal and maternal.

## Output:
## 1. Chromosome, Allele position, allele1 in Male, allele2 in Male (being '.' when male is homo), allele1 from Female, allele2 from Female (being '.' when female is homo).
## 2. 3 columns per sample:
##    a. allele inherited from male
##    b. allele inherited from female
##    c. total reads coverage of that position


usage = "\n\n%prog"
parser = OptionParser(usage=usage)

parser.add_option("-a", "--allele", dest="a", help="Input coverage filtered allele file.")
parser.add_option("-o", "--outputfile", dest="o", help="Output file.")

(opt, args) = parser.parse_args()

afile = open(opt.a) #open input coverage filtered call allele file
o = open(opt.o,'w')

## Write the header of now synthesized file ##
indices = []
totalcount = 0

line = afile.readline()
l = line.split('\t')
if line[0] == 'C':
    count = 0
    header = "Chrom\tPos\tM1\tM2\tF1\tF2"
    for m in l:
        if m.split('-')[0] == 'TotalCov':
            indices.append(count)
            lib = m[9:]
            header += '\tMale-'+lib+'\tFemale-'+lib+'\tTotalCov-'+lib
        count += 1
    o.write(header+'\n')

## Defining the inheritance of each SNP, and write the output file. 
for line in afile:
    l = line.split('\t')
    chrom = l[0]
    pos = l[1]
    m1 = l[2]
    m2 = l[3]
    f1 = l[4]
    f2 = l[5]
    text = chrom+'\t'+str(pos)+'\t'+m1+'\t'+m2+'\t'+f1+'\t'+f2

    for i in indices:
        Cov = l[i]
        totalcount += 1
        snp1 = l[i-2]
        snp2 = l[i-1]
        if snp2 == '.':
            text += '\t'+snp1+'\t'+snp1+'\t'+str(Cov)
        elif snp2 != '.':
            if m2 == '.' and m1 == snp1:
                text += '\t'+snp1+'\t'+snp2+'\t'+str(Cov)
            elif m2 == '.' and m1 == snp2:
                text += '\t'+snp2+'\t'+snp1+'\t'+str(Cov)
            elif f2 == '.' and f1 == snp1:
                text += '\t'+snp2+'\t'+snp1+'\t'+str(Cov)
            elif f2 == '.' and f1 == snp2:
                text += '\t'+snp1+'\t'+snp2+'\t'+str(Cov)
            elif m2 != '.' and f2 != '.':
                maletb=[m1,m2]
                femaletb=[f1,f2]
                if ((snp1 in maletb) and (snp1 not in femaletb)):
                    text += '\t'+snp1+'\t'+snp2+'\t'+str(Cov)
                elif ((snp1 in femaletb) and (snp1 not in maletb)):
                    text += '\t'+snp2+'\t'+snp1+'\t'+str(Cov)
                elif ((snp2 in maletb) and (snp2 not in femaletb)):
                    text += '\t'+snp2+'\t'+snp1+'\t'+str(Cov)
                elif ((snp2 in femaletb) and (snp2 not in maletb)):
                    text += '\t'+snp1+'\t'+snp2+'\t'+str(Cov)
                else:
                    text += '\t.'*2+'\t'+str(Cov)
            else:
                text += '\t.'*2+'\t'+str(Cov)
    text += '\n'
    o.write(text)

o.close()
afile.close()
                    
                


            


        

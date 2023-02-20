#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 18:10:43 2019
@author: wendy
"""

from optparse import OptionParser

#male is SO3615L
#female1 is SO5465L
#female2 is SO5985L

usage = ""
parser = OptionParser(usage=usage)
parser.add_option("-f", "--input_file", dest="f", help="input parsed mpileup file.")

(opt, args) = parser.parse_args()

mpufile = opt.f

def hetero_justify(SNP1,SNP2,SNP3,COV):
    if COV != '.' and int(COV) >= 20 :
        if SNP2 == '.':
            first = SNP1.split('_')
            nuleo = [first[0]]
            return nuleo
        elif SNP3 != '.':
            nuleo = ['NA']
            return nuleo
        elif SNP2 != '.':
            first = SNP1.split('_')
            second = SNP2.split('_')
            if float(first[1]) < 25:
                nucleo = [second[0]]
                return nucleo
            elif float(first[1]) > 75:
                nucleo = [first[0]]
                return nucleo
            elif float(first[1]) <= 75 and float(first[1]) >= 25:
                nucleo = [first[0],second[0]]
                return nucleo
        elif SNP1 == '.':
            nucleo = ['NA']
            return nucleo
    else:
        nucleo = ['NA']
        return nucleo

#list1 is male(homo)+female1(heter)
#list2 is male(homo)+female2(heter)
#list3 is male(heter)+female1(not same heter)
#list4 is male(heter)+female2(not same heter)
list1=[]
list2=[]
list3=[]
list4=[]

#open input mpileup file
mpu = open(mpufile,'r')

while 1:
    l = mpu.readline()
    if l == '':
        break
    ls = l.split('\t')
    if ls[0] == 'Chrom':
        continue
    else:
        covm = ls[14]
        covm1= covm.split('\n')
        nuc_m = hetero_justify(ls[11],ls[12],ls[13],covm1[0])
        if len(nuc_m) == 1 and nuc_m[0] != 'NA':
            nuc_f1 = hetero_justify(ls[3],ls[4],ls[5],ls[6])
            nuc_f2 = hetero_justify(ls[7],ls[8],ls[9],ls[10])
            #list1 and list2 will have overlapped SNP when female1 and female2 both are hetero and male is homo
            if len(nuc_f1) == 2:
                inls = [ls[0],ls[1],'m' + '_' + nuc_m[0], covm1[0], 'f' + '_' + nuc_f1[0] + '/' + nuc_f1[1], ls[6]]
                list1.append(inls)
            if len(nuc_f2) == 2:
                inls = [ls[0],ls[1],'m' + '_' + nuc_m[0],covm1[0], 'f' + '_' + nuc_f2[0] + '/' + nuc_f2[1], ls[10]]
                list2.append(inls)
        elif len(nuc_m) == 2:
            nuc_f1 = hetero_justify(ls[3],ls[4],ls[5],ls[6])
            nuc_f2 = hetero_justify(ls[7],ls[8],ls[9],ls[10])
            if len(nuc_f1) == 1 and nuc_f1[0] != 'NA':
                inls = [ls[0],ls[1],'m' + '_' + nuc_m[0] + '/' + nuc_m[1],covm1[0],'f' + '_' + nuc_f1[0],ls[6]]
                list3.append(inls)
            elif len(nuc_f1) == 2:
                if (nuc_f1[0] == nuc_m[0] and nuc_f1[1] == nuc_m[1]) or (nuc_f1[0] == nuc_m[1] and nuc_f1[1] == nuc_m[0]):
                    continue
                else:
                    inls = [ls[0],ls[1],'m' + '_' + nuc_m[0] + '/' + nuc_m[1],covm1[0],'f' + '_' + nuc_f1[0] + '/' + nuc_f1[1],ls[6]]
                    list3.append(inls)      
            if len(nuc_f2) == 1 and nuc_f2[0] != 'NA':
                inls = [ls[0],ls[1],'m' + '_' + nuc_m[0] + '/' + nuc_m[1],covm1[0],'f' + '_' + nuc_f2[0],ls[10]]
                list4.append(inls)
            elif len(nuc_f2) == 2:
                if (nuc_f2[0] == nuc_m[0] and nuc_f2[1] == nuc_m[1]) or (nuc_f2[0] == nuc_m[1] and nuc_f2[1] == nuc_m[0]):
                    continue
                else:
                    inls = [ls[0],ls[1],'m' + '_' + nuc_m[0] + '/' + nuc_m[1],covm1[0],'f' + '_' + nuc_f2[0] + '/' + nuc_f2[1],ls[10]]
                    list4.append(inls)

mpu.close()

o1 = open('SNP-list1.txt','w')
o2 = open('SNP-list2.txt','w')
o3 = open('SNP-list3.txt','w')
o4 = open('SNP-list4.txt','w')
for item in list1:
    o1.write(item[0] + '\t' + item[1] + '\t' + item[2] + '\t' + item[3] + '\t' + item[4] + '\t' + item[5] + '\n')
o1.close()

for item in list2:
    o2.write(item[0] + '\t' + item[1] + '\t' + item[2] + '\t' + item[3] + '\t' + item[4] + '\t' + item[5] + '\n')
o2.close()         
            
for item in list3:
    o3.write(item[0] + '\t' + item[1] + '\t' + item[2] + '\t' + item[3] + '\t' + item[4] + '\t' + item[5] + '\n')
o3.close()
            
for item in list4:
    o4.write(item[0] + '\t' + item[1] + '\t' + item[2] + '\t' + item[3] + '\t' + item[4] + '\t' + item[5] + '\n')
o4.close()
                


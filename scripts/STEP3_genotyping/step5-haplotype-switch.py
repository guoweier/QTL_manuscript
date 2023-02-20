from __future__ import division
import os, sys, math, time
from optparse import OptionParser
from collections import OrderedDict

#this script is used to correct the haplotype based on the potential switch points from haplotype-varify.py script. 

#input: a. list from haplotype-varify.py output. The potential switch point of previous haplotype.
#       b. the initial haplotype
#output: switched haplotype list



usage = "\n\n%prog"
parser = OptionParser(usage=usage) 

parser.add_option("-H", "--hapotype", dest="H", help="input initial non-switched haplotype list")
parser.add_option("-s", "--switch", dest="s", help="input potential switch point list from haplotype varify script")
parser.add_option("-o", "--outputfile", dest="o", help="output file")

(opt, args) = parser.parse_args()

hfile = open(opt.H)
sfile = open(opt.s)
o = open(opt.o, "w")

shead = sfile.readline()
hhead = hfile.readline()
switchdic = {}
haplodic = {}

for line in sfile:
	line = line[:-2]
	ln = line.split("\t")
	chrom = int(ln[0][3:])
	if chrom not in switchdic:
		switchdic[chrom] = []
	switchdic[chrom] += [ln[1]]
print(switchdic)

for line in hfile:
	line = line[:-1]
	ln = line.split("\t")
	chrom = int(ln[0][3:])
	if chrom not in haplodic:
		haplodic[chrom] = []
	haplodic[chrom] += [[ln[1],ln[2],ln[3]]]

for c in switchdic:
	for i in range(len(switchdic[c])):
		if (i % 2) == 0:
			j = i+1
			pos1 = switchdic[c][i]
			if j == len(switchdic[c]):
				for k in range(len(haplodic[c])):
					if int(haplodic[c][k][0]) > int(pos1):
						m1 = haplodic[c][k][1]
						m2 = haplodic[c][k][2]
						haplodic[c][k][1] = m2
						haplodic[c][k][2] = m1							
			elif j < len(switchdic[c]):
				pos2 = switchdic[c][j]
				for k in range(len(haplodic[c])):
					if int(haplodic[c][k][0]) > int(pos1) and int(haplodic[c][k][0]) <= int(pos2):
						m1 = haplodic[c][k][1]
						m2 = haplodic[c][k][2]
						haplodic[c][k][1] = m2
						haplodic[c][k][2] = m1

#output new switched list
o.write(hhead)
for item in haplodic:
	if len(str(item)) == 1:
		chrm = "Chr0"+str(item)
	elif len(str(item)) == 2:
		chrm = "Chr"+str(item)
	for i in range(len(haplodic[item])):
		pos = haplodic[item][i][0]
		m1 = haplodic[item][i][1]
		m2 = haplodic[item][i][2]
		text = chrm+"\t"+pos+"\t"+m1+"\t"+m2
		o.write(text+"\n")

hfile.close()
sfile.close()
o.close()













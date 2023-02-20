# Weier Guo python
# Created: 06/10/2021
# Updated: 05/10/2022 (record bin positions to be start and end, and exhibit in: Chr_Start_End)

import os, sys, math, time
from optparse import OptionParser
from collections import OrderedDict

#This script convert bin-by-SNP# output file into read.cross csv file (genotype) to be ready for R/qtl package.


usage = "\n\n%prog"
parser = OptionParser(usage=usage)

parser.add_option("-f", "--file", dest="f", help="Input bin-by-SNP# output file")
parser.add_option("-o", "--outputfile", dest="o", help="Output file.")

(opt, args) = parser.parse_args()

#input parameters
f = open(opt.f)
o = open(opt.o, "w")

Libs = {}
Markers = {}
count = 0

#read input file header
#record sample names and their haplotype info columns
line = f.readline()
line = line.split("\n")[0]
ln = line.split("\t")
for i in ln:
	if i.split("-")[0] == "Haplotype":
		lib = i[10:]
		Libs[lib,count] = {}
	count += 1

for line in f:
	line = line.split("\n")[0]
	ln = line.split("\t")
	chrom = ln[0]
	start = ln[1]
	end = ln[2]
	startpos = ln[3]
	endpos = ln[4]
	marker = chrom.split("Chr")[1]+"_"+start+"_"+end
	if len(start) < 7:
		markerS = '0'*(7-len(start))+start
	else:
		markerS = start
	markerNum = int(chrom.split("Chr")[1]+str(markerS))
	Markers[markerNum] = [chrom,start,end,marker,startpos,endpos]
	for i in Libs:
		if ln[i[1]] == "1":
			Libs[i][markerNum] = "A"
		elif ln[i[1]] == "2":
			Libs[i][markerNum] = "B"
		elif ln[i[1]] == ".":
			Libs[i][markerNum] = "NA"

f.close()

#sort dictionary
sortedmarkers = Markers.keys()
sortedmarkers.sort()
sortedlibs = Libs.keys()
sortedlibs.sort()

#set output header
header = "Chrom,Start,End,Marker,Pos"
for Sample in sortedlibs:
	header += ","+Sample[0]
o.write(header+"\n")

#write data into output file
for i in sortedmarkers:
    Chrom = str(Markers[i][0])
    Start = str(Markers[i][1])
    End = str(Markers[i][2])
    Marker = str(Markers[i][3])
    Pos = Chrom+'_'+str(Markers[i][4])+'_'+str(Markers[i][5])
    text = Chrom+','+Start+','+End+','+Marker+','+Pos

    for Sample in sortedlibs:
    	haplo = Libs[Sample][i]
    	text += ","+str(haplo)
    o.write(text+'\n')

o.close()








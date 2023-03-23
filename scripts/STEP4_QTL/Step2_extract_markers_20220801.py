# Weier Guo Python
# 08/01/2022 create
# 11/18/2022 update: extract more valuable information from the new combined markers. To do this, I get alleles based on the frequency of NA in combined alleles.
# 12/05/2022 update: remove the 11/18/2022 update 

import os, sys, math
from optparse import OptionParser
import itertools
from collections import OrderedDict
import statistics

usage = ""
parser = OptionParser(usage=usage)
parser.add_option("-s", "--standard_file", dest="s", help="Input standard marker file.")
parser.add_option("-g", "--input_file", dest="g", help="Input file with alleles information covering the genome.")
parser.add_option("-o", "--output", dest="o", help="Output file.")

(opt, args) = parser.parse_args()

### PARAMETERS DEFINE ###
s = open(opt.s,"r")
g = open(opt.g,"r")
o = open(opt.o,"w")
New_markers = {}


# read standard marker file head
s_head = s.readline()
s_head = s_head.split("\n")[0]
s_hd = s_head.split("\t")[1:]

# read genome marker file head
g_head = g.readline()
g_head = g_head.split("\n")[0]
g_hd = g_head.split("\t")[1:]

# function: count each variable numbers in the list
def var_num(L):
	result = {}
	L_un = list(set(L))
	for x in L_un:
		count = 0
		for i in L:
			if i == x:
				count += 1
		result[x] = count
	return result


### SET THE HEADER OF NEW FILE ###
# get each standard marker info
for i in range(len(s_hd)):
	smarker = s_hd[i]
	s_chr = smarker.split("_")[0]
	s_start = int(smarker.split("_")[1])
	s_end = int(smarker.split("_")[2])
	New_markers[i] = []
	
	# search genome marker position of standard marker
	for j in range(len(g_hd)):
		gmarker = g_hd[j]
		g_chr = gmarker.split("_")[0]
		g_start = int(gmarker.split("_")[1])
		g_end = int(gmarker.split("_")[2])
		if g_chr == s_chr:
			if s_start >= g_start and s_start <= g_end:
				New_markers[i] += [j]		# add marker start position
			if s_end >= g_start and s_end <= g_end:
				New_markers[i] += [j]		# add marker end position
				break 		# stop the current loop
				
	# arrange marker positions
	if New_markers[i][0] == New_markers[i][1]:
		New_markers[i] = [New_markers[i][0]]			# if start and end in the same genome marker, get only one marker position
	elif New_markers[i][1] - New_markers[i][0] > 1:
		New_markers[i] = [num for num in range(New_markers[i][0],New_markers[i][1]+1)]		# if start and end not in adjacent genome marker, get all the marker positions between them

New_markers = OrderedDict(sorted(New_markers.items()))

# write the header of output file
ohd = "Old.Name"
for marker in s_hd:		# the output file is using standard markers
	ohd += "\t"+marker 
o.write(ohd+"\n")

### GET ALLELES IN GENOME MARKER FILE ###
# read alleles info from genome marker file
for line in g:
	line = line.split("\n")[0]
	sample = line.split("\t")[0]
	ln = line.split("\t")[1:]
	text = sample 
	for key in New_markers:
		alleles = []
		for i in New_markers[key]:
			alleles.append(ln[i])
		alleles_count = var_num(alleles)
		if "NA" in alleles:
			allele = "NA"
		else:
			Alleles = list(map(lambda x: round(float(x),1), alleles))
			allele = round(statistics.mean(Alleles),1)
		text += "\t" + str(allele)
	o.write(text+"\n")

### END ###
g.close()
s.close()
o.close()







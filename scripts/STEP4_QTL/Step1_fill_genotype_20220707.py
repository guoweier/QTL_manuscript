# Weier Guo Python
# 07/07/2022 created
# 09/17/2022 updated: remove alleles which belongs to the skipped add-on markers because of the adjacency of two original markers. 
# 11/15/2022 updated: F1 genotype does not have marker numbers row. Remove the first header reading for F1 genotypes. 


import os, sys, math
from optparse import OptionParser
import itertools
from collections import OrderedDict

usage = ""
parser = OptionParser(usage=usage)
parser.add_option("-f", "--input_file", dest="f", help="Genetic map file with markers physical position.")
parser.add_option("-g", "--genome_size", dest="g", help="Input file of each chromosome size.")
parser.add_option("-o", "--output", dest="o", help="Output file.")

(opt, args) = parser.parse_args()

#### FUNCTION AND PARAMETER DEFINE ####
# define parameters
f = open(opt.f,"r") 
g = open(opt.g,"r") 
o = open(opt.o,"w")
Genome_mk = {}				# dictionary for all markers
Genomes = OrderedDict()		# dictionary for each chromosome size
Ori_mk = OrderedDict()		# dictionary for original markers
Mk_num = {}					# dictionary for marker number with corresponded chromosome

# define pairwise() function
def pairwise(iterable):
	a,b = itertools.tee(iterable)
	next(b, None)
	return zip(a,b)

# concatenate two lists by alternating elements
def knit_lists(a,b):
	c = itertools.zip_longest(a,b,fillvalue='-')
	result = []
	for value in c:
		result.append(value[0])
		result.append(value[1])
	return result

# get alleles for newly synthesized marker
def get_allele(a,b):
	if a == b:
		new = str(a) 
	elif a == "NA" and b != "NA":
		new = b 
	elif a != "NA" and b == "NA":
		new = a 
	else:
		new = "NA"
	return new


#### SET HEADER WITH MARKERS COVER THE GENOME ####
# adjust each chromosome size
for line in g:
	if line[0] == "C":
		line = line.split("\n")[0]
		ln = line.split("\t")
		c = int(ln[0][3:])
		Genomes[c] = ln[1]

# get markers physical position list
# mname_head = f.readline()
mpos_head = f.readline()
mpos_head = mpos_head.split("\n")[0]
mpos_hd = mpos_head.split(",")[1:]

# separate original markers by chromosome
for marker in mpos_hd:
	chrom = int(marker.split("_")[0][3:])
	if chrom not in Ori_mk:
		Ori_mk[chrom] = []
		index = mpos_hd.index(marker)
		Mk_num[index] = chrom
	Ori_mk[chrom] += [marker]


# get adjacent new markers
mpos_pair = pairwise(mpos_hd)	# get paired markers
for item in mpos_pair:
	m1 = item[0]				# get marker1
	m2 = item[1]				# get marker2
	chr1 = m1.split("_")[0]		# get chr for marker1
	chr2 = m2.split("_")[0]		# get chr for marker2
	start1 = m1.split("_")[1]	# get start pos for marker1
	start2 = m2.split("_")[1]	# get start pos for marker2
	end1 = m1.split("_")[2]		# get end pos for marker1
	end2 = m2.split("_")[2]		# get end pos for marker2

	if chr1 != chr2:
		continue
	else:
		if int(chr1[3:]) not in Genome_mk:
			Genome_mk[int(chr1[3:])] = []
			str_marker = chr1+"_"+"1"+"_"+str(int(start1)-1)
			Genome_mk[int(chr1[3:])] += [str_marker]
			for marker in Ori_mk[int(chr1[3:])]:
				Genome_mk[int(chr1[3:])] += [marker]
		# add newly snythesized marker
		adj_marker = chr1+"_"+str(int(end1)+1)+"_"+str(int(start2)-1)
		Genome_mk[int(chr1[3:])] += [adj_marker]

# add the last marker for each chromosome
for key in Genome_mk:
	emark_chr = Genome_mk[key][-1].split("_")[0]
	emark_str = str(int(Ori_mk[key][-1].split("_")[2])+1)
	emark_end = Genomes[key]
	end_marker = emark_chr+"_"+emark_str+"_"+emark_end
	Genome_mk[key] += [end_marker]


# sort markers in each chromosome
Genome_mk = OrderedDict(sorted(Genome_mk.items()))
for key in Genome_mk:
	Genome_mk[key].sort(key = lambda m: int(m.split("_")[1]))

# write the title for output file
ohead = "Old.Name"
for key in Genome_mk:
	for marker in Genome_mk[key]:
		ohead += "\t"+marker
o.write(ohead+"\n")

#### GET SAMPLES INFORMATION ####
# get alleles for newly synthesized markers for samples
for line in f:
	line = line.split("\n")[0]
	ln = line.split(",")[1:]
	Sample = line.split(",")[0]

	ori_al = {}			# dictionary for original sample alleles
	adj_al = {} 		# dictionary for newly synthesized alleles
	genome_al = {} 		# dictionary for concatenating original and synthesized alleles

	# separate original alleles by chromosome
	for i in range(len(ln)):
		if i in Mk_num:
			chrom = Mk_num[i]
			ori_al[chrom] = [ln[i]]
		else:
			ori_al[chrom] += [ln[i]]

	ori_al = OrderedDict(sorted(ori_al.items()))

	# get adjacent new alleles
	for chrom in ori_al:
		adj_ini = list(pairwise(ori_al[chrom]))
		for pair in adj_ini:
			if pair[1] != None:
				if chrom not in adj_al:
					adj_al[chrom] = []
				adj_al[chrom] += [get_allele(pair[0],pair[1])]
		adj_al[chrom] += [adj_ini[-1][1]]

	# knit original and newly synthezied allele lists
	for chrom in ori_al:
		genome_al[chrom] = [ori_al[chrom][0]]
		comb_al = knit_lists(ori_al[chrom],adj_al[chrom])
		for item in comb_al:
			genome_al[chrom] += [item]

	genome_al = OrderedDict(sorted(genome_al.items()))


	# write alleles into text
	text = Sample
	for chrom in genome_al:
		for allele in genome_al[chrom]:
			text += "\t"+str(allele)
	o.write(text+"\n")


f.close()
g.close()
o.close()






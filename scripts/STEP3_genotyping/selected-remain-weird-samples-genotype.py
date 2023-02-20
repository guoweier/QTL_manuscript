import os, sys, math, time
from optparse import OptionParser

usage = "\n\n%prog"
parser = OptionParser(usage=usage)

parser.add_option("-s", "--selected", dest="s", help="Input selected het haplotype file")
parser.add_option("-o", "--outputfile", dest="o", help="Output file.")

(opt, args) = parser.parse_args()

#open files#
f = open(opt.s)
o = open(opt.o, "w")

indices = []
count = 0

fhead = f.readline()
fheader = fhead.split("\t")
delsam = ["4","19","27","33","47","49","53","69","90","92","95","113","126","150","172","182","193","199","240","255","277","304","308","320","345","369","392","428","433","443","445","450","477","479","483","488","32","58","60","71","80","118","170"]
for item in fheader:
    if item.split("-")[0] == "GWR" or item.split("-")[0] == "XXX":
        if item.split("-")[2] in delsam:
            indices.append(count)
    count += 1

#output file header#
header = "Chrom\tPos\tM1\tM2\tF1\tF2"
for a in indices:
    factor = fheader[a]
    header += "\t"+factor
o.write(header+"\n")

for line in f:
    ln = line.split("\t")
    wline = ln[0]+"\t"+ln[1]+"\t"+ln[2]+"\t"+ln[3]+"\t"+ln[4]+"\t"+ln[5]
    for a in indices:
        factor = ln[a]
        wline += "\t"+factor
    o.write(wline+"\n")

f.close()
o.close()

#Usage: python pileup2seq.py -f AA_1013filmultiN.dat -p bwtconsensus.pileup
#Function: Extract 101 bp flanking region (50 bp upstream and downstream) of a SNP using consensus pileup file. The output file is a CSV format file designated for non-model organism in ADT submission.
#This script is written by Melissa Wong (melissawongukm@gmail.com)

from optparse import OptionParser
from os.path import splitext
import re
p = OptionParser()
p.add_option("-f", dest="f1",help="SNPcount.out file")
p.add_option( "-p", dest="p1", help="consensus pileup file")
(opt,args) = p.parse_args()
f1 = open(opt.f1)
p1 = open(opt.p1)
out = open(splitext(opt.f1)[0]+'probe.csv','w')
out.write("Locus_Name,Target_Type,Sequence,Chromosome,Coordinate,Genome_Build_Version,Source,Source_Version,Sequence_Orientation,Plus_Minus\n")
idict = {'A/C/T/G':'X','A/C/G':'V','C/T/G':'B','A/C/T':'H','D':'A/T/G',\
         'T/G':'K','C/G':'S','A/T':'W','A/C':'M','C/T':'Y','A/G':'R'}
tagset = set()
rlist = []
allSets = {}
for line in f1:
    line = line.strip()
    tags = line.split('\t')
    tagset.add(tags[0])
    columns = tags[:2]
    position = int(tags[1])
    snplist = "\t".join(columns)
    allSets[snplist] = set()   
    allSets[snplist] = allSets[snplist].union(range(position-50, position+51))     
f1.close()
for r in p1:
    s = r.split()
    tag = r.split()[0]
    if tag in tagset:
        rlist.append(s)
p1.close()
klist = allSets.keys()
for key in klist:
    key1 = key.split()[0]
    key2 = key.split()[1]
    bases = []
    for r in rlist:
	if r[0] == key1 and int(r[1]) in allSets[key]:
 	    if r[1]== key2:
  	        bases.append('[')            
  	        bases.append(r[2])
  	        bases.append(']') 
	    elif r[2] in idict:
                bases.append(idict[r[2]])
  	    else: 
	   	bases.append(r[2])
    allbases = "".join(bases)
    out.write("%s_%s,SNP,%s,0,0,0,unknown,0,unknown,Plus\n" % (key1, key2, allbases)) 
out.close()

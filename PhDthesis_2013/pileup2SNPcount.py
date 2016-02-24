#Usage: python pileup2SNPcount.py -f AArefseq_bwtAA6_fil50bs.dat
#Function: Calculate the number of non-reference allele, number of third allele, minor allele frequency and major allele frequency from Samtools pileup file. This script is modified to score missed SNPs that do not share the same allele with reference sequences. Minor allele frequency and major allele frequency are put in the correct order after sorting in list.
#This script is written by Melissa Wong (melissawongukm@gmail.com)

from __future__ import division
from optparse import OptionParser
from os.path import splitext
import re

p = OptionParser()
p.add_option("-f", dest="f1",help="file")
(opt,args) = p.parse_args()
f1 = open(opt.f1)
out = open(splitext(opt.f1)[0]+'_SNPselect.out','w')

ref = re.compile (r"([.,])",re.IGNORECASE)
l = re.compile (r"([.,ACGT])",re.IGNORECASE)
IUB = {'K':re.compile (r"([GT])",re.IGNORECASE),'S':re.compile (r"([GC])",re.IGNORECASE),\
       'W':re.compile (r"([AT])",re.IGNORECASE),'M':re.compile (r"([AC])",re.IGNORECASE),\
       'Y':re.compile (r"([CT])",re.IGNORECASE),'R':re.compile (r"([AG])",re.IGNORECASE)}
idict = {'K':['T','G'],'S':['C','G'],'W':['A','T'],'M':['A','C'],'Y':['C','T'],'R':['A','G']}
alel = {'C':re.compile (r"([Cc])"),'T':re.compile (r"([Tt])"),'A':re.compile (r"([Aa])"),'G':re.compile (r"([Gg])")}


for line in f1:
    column = line.split("\t")
    read = column[8]
    consensus = column[3]
    reflen = len(ref.findall(read))
    llen = len(l.findall(read)) 
    if reflen < 2:
        x = idict[consensus]
     	list1 = []
     	length = 0
	for i in x:
            ylen = len(alel[i].findall(read)) 
            freq = ylen/llen
            list1.append("%.2f" % freq)
            list1.append(ylen)
            length += ylen
        list1.sort()
        NR = list1[0]
	trialel = llen - length
	out.write("\t".join(column[:6]) + "\t")
	out.write("\t".join(list1[2:]))
	out.write("\t%d\t%d\t"%(trialel,NR))
	out.write("\t".join(column[7:]))
    else: 
	list2 = [] 
    	zlen = len(IUB[consensus].findall(read))
    	freq1 = reflen/llen
    	freq2 = zlen/llen    
    	trialel = llen - zlen - reflen
	list2.append("%.2f" % freq1)
	list2.append("%.2f" % freq2)
	list2.sort()
    	out.write("\t".join(column[:6]) + "\t")
	out.write("\t".join(list2))
    	out.write("\t%d\t%d\t"%(trialel, zlen))
    	out.write("\t".join(column[7:]))
out.close()

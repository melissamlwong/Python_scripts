#Usage: python pileup2consensus_bwt.py -f AArefseq_merge4bwt.pileup
#Function: Find the final consensus sequences from the consensus sequences of four Acacia parents. 
#This script is written by Melissa Wong (melissawongukm@gmail.com)

from optparse import OptionParser
from os.path import splitext
import re
p = OptionParser()
p.add_option("-f", dest="f1",help="Consensus sequences of four Acacia parents in pileup format")
(opt,args) = p.parse_args()
f1 = open(opt.f1)
out = open(splitext(opt.f1)[0]+'_consensus.pileup','w')
flag = False
consensus = []
IUPAC = ['K','W','S','M','Y','R']
IUB = set(IUPAC)
idict = {'X':['A','C','T','G'],'V':['A','C','G'],'B':['C','T','G'],'H':['A','C','T'],'D':['A','T','G'],\
         'K':['T','G'],'S':['C','G'],'W':['A','T'],'M':['A','C'],'Y':['C','T'],'R':['A','G']}
#use X to distinguish difference with N
def find_key(dic, val):
   return [k for k, v in idict.iteritems() if v == val][0]
for line in f1:
    l = line.strip()
    columns = l.split("\t")
    reads = columns[3:]
    pairAA = columns[3:5]
    pairAM = columns[5:]
    total = set(reads)
    total.discard('N')
    snp = total.intersection(IUB)
    totalist = list(total)
    total_str = "".join(reads)
    if reads.count('N') >= 3:
        consensus = 'N'
    elif pairAA.count('N') == 2 or pairAM.count('N') == 2:
        consensus = 'N'
    elif len(total) == 1:
 	consensus = list(total)
    elif len(snp) > 0:
        for i in IUPAC:
	    if i in snp:
                total|=set(idict[i])
                total.remove(i)
                consensus = list(total)
    else:
    	I = find_key(idict, totalist)
        consensus = I
    out.write("\t".join(columns) + "\t" + "/".join(consensus) + "\n")
out.close()

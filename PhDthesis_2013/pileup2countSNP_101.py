#Usage: pileup2countSNP_101.py -f AA_1013filmultiN.dat -p bwtconsensus.pileup
#Function: count the number of secondary SNPs, triallelic SNP, interspecific SNP and Ns within 101 bp flanking region of a SNP using consensus pileup file. The performance of this script has been improved. 
#This script is written by Melissa Wong (melissawongukm@gmail.com)

from optparse import OptionParser
from os.path import splitext
import re
p = OptionParser()
p.add_option("-f", dest="f1",help="Filtered SNPs pileup file")
p.add_option( "-p", dest="p1", help="consensus pileup file")
(opt,args) = p.parse_args()
f1 = open(opt.f1)
p1 = open(opt.p1)
out = open(splitext(opt.f1)[0]+'count.out','w')
out.write("contig\tposition\tlength\tsnp\ttsnp\tIUB\tN\n")
idict = {'isnp':re.compile("[KWSMYRXVBHD]"), 'snp':re.compile("[ACGT]/[ACGT]"),\
         'tsnp':re.compile("/[ACGT]/"), 'N':re.compile("N")}
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
    allSets[snplist] = allSets[snplist].union(range(position-50, position+51)) #101
f1.close()
for r in p1:
    s = r.split()
    tag = s[0]
    if tag in tagset:
        rlist.append(s)
p1.close()
klist = allSets.keys()
for key in klist:
    key1 = key.split()[0]
    length = 0
    bases = []
    for r in rlist:
	if r[0] == key1 and int(r[1]) in allSets[key]:
	    length += 1 #total length 101
  	    bases.append(r[2])
    allbases = "".join(bases)
    a = len(idict['snp'].findall(allbases))
    b = len(idict['tsnp'].findall(allbases))
    c = len(idict['isnp'].findall(allbases))
    e = len(idict['N'].findall(allbases))
    f = a - b #true snps after removing triallelic snp; Assuming no quad-allelic snps
    out.write("%s\t%d\t%d\t%d\t%d\t%d"%(key,length,f,b,c,e)+ "\n")
out.close()

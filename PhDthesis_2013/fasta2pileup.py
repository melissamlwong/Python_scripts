#Usage: python fasta2pileup.py -f AArefseq_merge.fa
#Function: Convert fasta file to pileup format. This script is written by Zhang Di for Melissa Wong (melissawongukm@gmail.com). 

from optparse import OptionParser
from os.path import splitext
p = OptionParser()
p.add_option("-f", dest="fa", help="fasta file.")
(opt,args) = p.parse_args()
fa = opt.fa
out = open(splitext(opt.fa)[0]+'_pileup','w')
for line in open(fa):
    if line.startswith('>'):
        tag = line.strip()[1:]
        index = 1
    else:
        for bp in line.strip():
            out.write("%s\t%d\t%s\n"%(tag,index,bp))
            index += 1
out.close()

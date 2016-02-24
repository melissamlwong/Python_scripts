#!/usr/bin/env python

from Bio import SeqIO
from os.path import splitext
from optparse import OptionParser
import sys, re
if len(sys.argv) < 2:
    sys.exit("\nUsage: python %s [filename.fa]\nCalculate number of fragments (default range: 5,000-15,000 bp) for single digestion\n" % sys.argv[0])
p = OptionParser() 
p.add_option("-f", action="store", type="string", dest="fasta", help="fasta file")
p.add_option("-e", action="store", type="string", dest="enzyme", help="single digestion enzyme sequence") 
p.add_option("-c", action="store", type="int", dest="adjust", default=int("0"),help="enzyme cut position")
p.add_option("-l", action="store", type="int", dest="min", default=int("5000"),help="lower limit of fragment size")
p.add_option("-u", action="store", type="int", dest="max", default=int("15000"),help="upper limit of fragment size")
p.add_option("-n", action="store", type="string", dest="name", default="unknown enzyme",help="enzyme name") #make this optional
(opt,args) = p.parse_args()
def calc_frag(a,b,c,d,e,f):
    t1=0;t2=0;t3=0
    for record in SeqIO.parse(f, "fasta"):
        tag = record.id
        sequence = str(record.seq)
        length = len(sequence)
        clist = [c.start()+b for c in re.finditer(str(c), sequence)] 
        clist.insert(0, int("0"))
        clist.append(length) 
        flist = [clist[i]-clist[i-1] for i in range(len(clist)) if i is not 0]
        s1 = [x for x in flist if x <d]
        s2 = [x for x in flist if x >=d and x<=e]
        s3 = [x for x in flist if x >e]
        t1+=len(s1);t2+=len(s2);t3+=len(s3)
    return "Total %s fragments in <%s,%s-%s,>%s kb is ...\t%s\t%s\t%s" % (a,d,d,e,e,t1,t2,t3)
if __name__ == '__main__':
    f1 = open(opt.fasta)
    print calc_frag(opt.name,opt.adjust,opt.enzyme,opt.min,opt.max,f1)
    f1.close()

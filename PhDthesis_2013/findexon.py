#Usage: python findexon.py -s Mtt3.5Acacia_2xmlfil.blast2 -t Mt3.5_IMGAGfilsortlen.gff2 -f 96plex_SNP.tagset
#Function: Predict the position of exon-intron boundary using blast xml alignments
##This script is written by Melissa Wong (melissawongukm@gmail.com)

from optparse import OptionParser
from os.path import splitext
parser = OptionParser()
parser.add_option( "-s", dest="fs", help="Parsed Blast results with XML alignments")
parser.add_option( "-t", dest="ft", help="Filtered Mt3.5 IMGAG GFF file")
parser.add_option( "-f", dest="f1", help="SNPs file")
(options, args) = parser.parse_args()
fs = open(options.fs)
ft = open(options.ft)
f1 = open(options.f1)
out = open(splitext(options.f1)[0]+'.exon2SNP','w')
out2 = open(splitext(options.f1)[0]+'.exondict','w')
exondict = {}; Acexondict = {}; tagset = set(); lendict = {}
def makelist(start, seq):
    list1 = []; list2 = ['A','C','T','G','N']
    y = int(start)
    for base in seq:
	if base in list2:
	   list1.append(str(y))
	   y += 1
	else:
	   list1.append('X') #Add "X" if there is a gap
    return list1
def counter(i, list1):
    while list1[i] == "X":
	i+=1
    return i
gff = ft.readlines()
for line in gff:
    Mtt = line.split()[0]; length = int(line.split()[3])
    lendict.setdefault(Mtt, []).extend([length])
    tagset.add(Mtt)
for Mtt in tagset:
    length = 1
    exon = []
    for item in lendict[Mtt]:
	length += item
	exon.append(length)
	length +=1
    exonboundary = [str(x) for x in exon[:-1]] + [str(x+1) for x in exon[:-1]]
    exonboundary.sort()
    exondict[Mtt] = exonboundary
    out2.write("%s\t%s\n" % (Mtt, exonboundary))
blast = fs.readlines()
for line in blast:
    Ac = line.split()[0]; Mtt = line.split()[1]; sseq = line.split()[13]; sstart = line.split()[8]; send = line.split()[9]
    if Mtt in exondict.keys():
        qlist = makelist(line.split()[6], line.split()[12])
        if int(sstart) < int(send):
            slist = makelist(sstart, sseq)
    	else:
	    slist = makelist(send, sseq[::-1])
	    slist.reverse()
        exon2pos = [i for i,x in enumerate(slist) if x in exondict[Mtt]]
	Ac2pos = [qlist[counter(i,qlist)] for i in exon2pos] #Find nearest position on the right if position = "X"
	Acexondict.setdefault(Ac, []).extend(Ac2pos)
SNP = f1.readlines()
for line in SNP:
    contig = line.split()[0]; snp = line.split()[1]
    if contig in Acexondict.keys():
        region = range((int(snp)-34),(int(snp)+34))
        position = [x for x in region if str(x) in Acexondict[contig]]
	if len(position) >0:
	    out.write("%s\t%s\n" % (contig, snp))
f1.close()
fs.close()
ft.close()
out.close()

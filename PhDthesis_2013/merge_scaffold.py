#Usage: python merge_scaffold.py -m AAAM_200fil.dat
#Function: Extract the longest orthologous sequences from two datasets using filtered Nucmer alignment results
#Automatically locate the sequence file using the destination from nucmer file
#Modified to extract only reference that match to one or more queries but not two or more references that match to one query. 
#This script is written by Zhang Di and modified by Melissa Wong (melissawongukm@gmail.com)

from optparse import OptionParser
from os.path import splitext
from collections import defaultdict
parser = OptionParser()
parser.add_option( "-m", dest="fm", help="the nucmer coords file.")
(options, args) = parser.parse_args()
fm = open(options.fm)
out = open(splitext(options.fm)[0]+'_merge.fa','w')
files = fm.readline().split()
(ref,query) = map(open,files)
matchdic = defaultdict(set)
refset = set()
for i in xrange(4):
    fm.readline()
for line in fm:
    segs = line.split('|');
    tags = segs[-1].split()
    matchdic[tags[1]].add(tags[0])
    refset.add(tags[0])
print 'total unique ref tags in coords',len(refset)
print 'total unique query tags in coords',len(matchdic)
fm.close()
for (k,v) in matchdic.items():
    if len(v)>=2:
        del matchdic[k]
print 'those unique match',len(matchdic)
flag = False
#for line in ref:
for line in query:
    if line.startswith('>'):
        flag = line[1:-1] in matchdic
    if flag:
        out.write(line)
ref.close()
query.close()
out.close()

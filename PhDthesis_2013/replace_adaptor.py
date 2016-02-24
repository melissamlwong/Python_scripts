from optparse import OptionParser
from os.path import splitext
from collections import defaultdict
parser = OptionParser()
parser.add_option( "-f", dest="fa", help="the fa file using as target.")
parser.add_option( "-m", dest="fm", help="the mummer output file.")

(options, args) = parser.parse_args()
fa = options.fa
fm = open(options.fm)

output = open(splitext(fa)[0]+'_Noadaptor.fa','w')
repdict = defaultdict(list)
for i in xrange(5):
    fm.readline()
for line in fm:
    segs = line.split('|');
    key = segs[-1].split()[-1]
    repdict[key].append(sorted(map(int,segs[1].split())))
print 'total unique tags in mummer out:',len(repdict)
print 'tags that have >1 replacements:',len([1 for v in repdict.values() if len(v)>1])
fm.close()
flag = False
for line in open(fa):
    if line.startswith('>'):
        tag = line[1:-1]
        flag = tag in repdict            
    elif flag:
        for (i,j) in repdict[tag]:        
            line=line[:i-1]+'N'*(j-i+1)+line[j:]
    output.write(line)
output.close()

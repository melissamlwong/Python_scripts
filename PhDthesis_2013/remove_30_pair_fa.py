from optparse import OptionParser
from os.path import splitext

p = OptionParser()
p.add_option("-x", dest="f1",help="the input fa file, p1")
p.add_option("-y", dest="f2",help="the input fa file, p2")
(opt,args) = p.parse_args()
f1 = open(opt.f1)
f2 = open(opt.f2)

out1s = open(splitext(opt.f1)[0]+'_30left1.fa','w')
out1b = open(splitext(opt.f1)[0]+'_30left.fa','w')
#out2s = open(splitext(opt.f2)[0]+'_17extractgood.fa','w')
out2b = open(splitext(opt.f2)[0]+'_30left.fa','w')

limit = 30
l1 = f1.readline()
l2 = f2.readline()
while (l1 != ''):
    if l1.startswith('>'):
        tag1 = l1
        tag2 = l2
    else:
        (i,j)=(max(map(len,l1.split('N')))<limit, max(map(len,l2.split('N')))<limit)
        if (i or j):
            if not i: out1s.write(tag1+l1)
            if not j: out1s.write(tag2+l2)            
        else:
            out1b.write(tag1+l1)
            out2b.write(tag2+l2)
    l1 = f1.readline()
    l2 = f2.readline()
out1s.close()
out1b.close()
out2b.close()
f1.close()
f2.close()


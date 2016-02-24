#This script was written by Melissa Wong (melissawongukm@gmail.com) to find reciprocal best hit from tab-delimited Blast results
#Output best hit for query
from collections import defaultdict
from os.path import splitext
import sys
if len(sys.argv) < 2:
    sys.exit('\nUsage: python %s [filename]\nReport reciprocal best hit from tab-delimited Blast results\n' \
    % sys.argv[0]) #python Reciprocalbesthit.py mt4genes_lenbacN.dat
out1 = open(splitext(sys.argv[1])[0]+'_pybesthit.dat','w')
out2 = open(splitext(sys.argv[1])[0]+'_rbesthit.dat','w')
f = open(sys.argv[1],"r")
fa = f.readlines()
f.close()
qdict = defaultdict(defaultdict); sdict = defaultdict(defaultdict); qbesthit = {}; sbesthit = {}
for line in fa: #Get total score for each query and each subject
    query,subject = line.split("\t")[0:2]
    score = int(round(float(line.strip().split("\t")[11])))
    if query not in qdict: #if no query 
        qdict[query][subject] = score
    elif query in qdict and subject not in qdict[query]: #if same query but not same subject
        qdict[query][subject] = score
    else: #if same query and subject add score
        total_score = qdict[query][subject] + score
        qdict[query][subject] = total_score
for query in qdict: #Now find best score
    maximum = 0
    for subject in qdict[query]:
        score = qdict[query][subject]
        if score > maximum:
            maximum = score
            qbesthit[query] = subject
    out1.write("%s\t%s\t%s\n" % (query, qbesthit[query],qdict[query][qbesthit[query]]))
for line in fa: #Repeat for reciprocal beshit
    query,subject = line.split("\t")[0:2]
    score = int(round(float(line.strip().split("\t")[11])))
    if subject not in sdict: #if no subject
        sdict[subject][query] = score
    elif subject in sdict and query not in sdict[subject]: #if same subject but not same query
        sdict[subject][query] = score
    else: #if same query and subject add score
        total_score = sdict[subject][query] + score
        sdict[subject][query] = total_score
for subject in sdict: #Now find best score
    maximum = 0
    for query in sdict[subject]:
        score = sdict[subject][query]
        if score > maximum:
            maximum = score
            sbesthit[subject] = query
#Find reciprocal besthit
for i in qbesthit:
    j = qbesthit[i]
    if sbesthit.has_key(j) and sbesthit[j] == i:
        out2.write("%s\t%s\n" % (i,j))
out1.close()
out2.close()

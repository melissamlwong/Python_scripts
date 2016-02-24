#Usage: python parsexml_biopython.py -f MtrcdsAc_blastxml.txt
#Function: parse particular information from blastn results in XML format
#This script is written by Melissa Wong (melissawongukm@gmail.com) based on Biopython documentation

from optparse import OptionParser
from os.path import splitext
from Bio.Blast import NCBIXML
p = OptionParser()
p.add_option("-f", dest="f1",help="blast results in xml format")
(opt,args) = p.parse_args()
result_handle = open(opt.f1)
out = open(splitext(opt.f1)[0]+'_2xml.blast','w')

blast_records = NCBIXML.parse(result_handle)
for record in blast_records:
    for alignment in record.alignments:
	Mtt = alignment.title 
    	for hsp in alignment.hsps:
	    out.write(("%s\t%s\t%d\t%s\t%s\t%s\t%s\n")%(record.query, Mtt.split(" ")[1], hsp.query_start, hsp.sbjct_start, hsp.bits, hsp.query, hsp.sbjct))
out.close()
result_handle.close()

#Predicts the longest Open Reading Frame (ORF) for each cDNA in the file. Return nucleotide and protein sequences of the ORF.
#This script was written by Melissa Wong (melissawongukm@gmail.com) using modified script from Biopython documentation

from Bio import SeqIO 
from Bio.Seq import reverse_complement, translate
from os.path import splitext
import sys

if len(sys.argv) < 2:
    sys.exit('\nUsage: python %s [filename]\nChoice of codon table and min protein length can be edited in the script\n' \
    % sys.argv[0])
out1 = open(splitext(sys.argv[1])[0]+'_cdspep.fa','w')
out2 = open(splitext(sys.argv[1])[0]+'_cdsnuc.fa','w')
table = 1 #standard table
min_pro_len = 1 #minimum protein length

def find_orfs_with_trans(seq, trans_table, min_protein_length):
    dictionary = {}
    seq_len = len(seq)
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start) #This is the ending index, by default its equal to the length of the string.
                if aa_end == -1: #Cannot find '*' anymore
                    aa_end = trans_len
                if aa_end-aa_start >= min_protein_length:  #Get start and end positions of all sequences regardless of strands
                    start = frame+aa_start*3
                    end = min(seq_len,frame+aa_end*3+3)
                    pro = trans[aa_start:aa_end]
                    gene = trans[start:end]
                    dictionary.update({len(pro) : [start, end, strand, pro, frame]}) #Hopefully, no max length is the same?
                aa_start = aa_end+1 #start looking for next * after the previous *
    return dictionary

for record in SeqIO.parse(sys.argv[1], "fasta"):
    tag = record.id
    dictionary = find_orfs_with_trans(record.seq, table, min_pro_len)
    key = sorted(dictionary.iterkeys())[-1] #find maximum protein length
    orf_list = dictionary[key]
    sequence = record.seq[orf_list[0]:orf_list[1]]
    rseq = record.seq.reverse_complement()
    rc = rseq[orf_list[0]:orf_list[1]]
    if orf_list[2] == 1:
        out1.write(">%s\n%s\n" % (tag, orf_list[3]))
        out2.write(">%s\n%s\n" % (tag, sequence))
    else:
        out1.write(">%s\n%s\n" % (tag, orf_list[3]))
        out2.write(">%s\n%s\n" % (tag, rc))
out1.close()
out2.close()

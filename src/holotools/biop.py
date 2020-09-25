#!/usr/bin/env python3
'''Just Make Biopython Work With Less Key Strokes'''
from Bio import SeqIO
def fdict(file):
    '''parse fasta file to dictionary'''
    d = {}
    with open(file,'r') as h:
        for r in SeqIO.parse(h,'fasta'):
            d[r.description]=str(r.seq)
    return d

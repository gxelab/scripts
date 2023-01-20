# -*- coding: utf-8 -*-
"""
shuffle RNA sequence in fasta format with ushuffle

parameters:
fasta_file: input fasta
prefix: prefix for output shuffled fasta
let: let = 2 means maintaining di-nucleotide frequency
n: shuffle how many times

example:
python fastaUshuffle.py test.fa res 2 10

use python 3.x only
------------------------------------------------
@author: zh (mt1022)
@date: Jan 14 2019
"""
import sys
from Bio import SeqIO

from ushuffle import shuffle, Shuffler

fasta_file, prefix, let, n = sys.argv[1:]
let = int(let)
n = int(n)

fasta = {}
for record in SeqIO.parse(fasta_file, 'fasta'):
    fasta[record.id] = record.upper()

outfh = {}
for i in range(n):
    outfile = prefix + (('_Ushuffle_%0' + str(len(str(n))) + 'd.fa') % i)
    fh = open(outfile, 'w')
    outfh[i] = fh

for record in fasta:
    seq = str(fasta[record].seq).encode('utf-8')
    seq_shuffler = Shuffler(seq, let)
    for i in range(n):
        new_seq = seq_shuffler.shuffle().decode('utf-8')
        print('>' + fasta[record].id, file=outfh[i])
        print(new_seq, file=outfh[i])

for i in range(n):
    outfh[i].close()


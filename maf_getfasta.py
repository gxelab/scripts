#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Get spliced multifasta from bed given bed intervals.

The results generated by the script give correct ref sequences
for 5' UTR regions of all the all the 71912 human transcripts
with 5' UTR.

requirement:
- biopython (>=1.72)

example usage:
    python maf_getfasta.py dir_of_maf hg38 input.bed output.mfa
"""

import sys
import csv
import re
import os
import traceback
from Bio.AlignIO import MafIO
from Bio import AlignIO
from Bio.SeqIO import FastaIO


def mglob(in_dir, pattern):
    """
    find files with a certain suffix in a dir
    """
    if not in_dir.endswith('/'):
        in_dir = in_dir + '/'
    res_files = [in_dir + i for i in os.listdir(in_dir) if i.endswith(pattern)]
    return res_files


def load_maf_index(maf_path, refsp):
    """
    load maf and mafidx files and reduction a dict{prefix:idx}
    """
    maf_files = mglob(maf_path, '.maf')
    maf_files.sort()
    # print(maf_files, file=sys.stderr)
    mafidx_files = mglob(maf_path, '.mafidx')
    mafidx_files.sort()
    # print(mafidx_files, file=sys.stderr)
    maf_prefix = [re.sub(r'.*(chr.*?)\.maf$', r'\1', i) for i in maf_files]
    print(maf_prefix, file=sys.stderr)
    mafidx_prefix = [re.sub(r'.*(chr.*?)\.mafidx$', r'\1', i) for i in mafidx_files]
    print(mafidx_prefix, file=sys.stderr)
    if set(maf_prefix) != set(mafidx_prefix):
        for i, j in zip(maf_files, maf_prefix):
            if j not in mafidx_prefix:
                print('....Error: No index found for maf: ' + i, file=sys.stderr)
        for i, j in zip(mafidx_files, mafidx_prefix):
            if j not in maf_prefix:
                print('....Error: No maf found for index: ' + i, file=sys.stderr)
        sys.exit(1)
    mafs = {i:j for i, j in zip(maf_prefix, maf_files)}
    mafidxes = {i:j for i, j in zip(mafidx_prefix, mafidx_files)}
    mafs_indexed = dict()
    for i in maf_prefix:
        mafs_indexed[i] = MafIO.MafIndex(mafidxes[i], mafs[i], '{}.{}'.format(refsp, i))
    return mafs_indexed


def load_splist(maf_path):
    """
    load ordered list of species prefix in maf
    """
    splist = []
    try:
        splist_path = '{}/{}'.format(maf_path, 'species_prefix_in_maf.txt')
        with open(splist_path, 'rt') as f:
            for line in f:
                splist.append(line.rstrip())
    except:
        print('Warning: each fasta block is not ordered due to the following error!', file=sys.stderr)
        print(traceback.format_exc(), file=sys.stderr)
    return splist


def parse_bed(mafs_indexed, splist, bed_path, out_path=None):
    """
    get spliced bed for each row of a bed file (12col or 6col)
    """
    if out_path is None:
        outh = sys.stdout
    else:
        outh = open(out_path, 'wt')
    with open(bed_path, 'rt') as bed:
        bed_reader = csv.reader(bed, delimiter='\t')
        for row in bed_reader:
            seq_name = row[3]
            chrom = row[0]
            if not chrom in mafs_indexed:
                print('{}: chromosome not found!'.format(chrom), file=sys.stderr)
                continue
            strand = 1 if row[5] is '+' else -1
            if len(row) >= 12:
                block_sizes = [int(i) for i in row[10].rstrip(',').split(',')]
                block_starts = [int(i) for i in row[11].rstrip(',').split(',')]
                region_starts = [i + int(row[1]) for i in block_starts]
                region_ends = [i + j for i, j in zip(region_starts, block_sizes)]
            else:  # not bed12 format, treat as 6 column bed
                region_starts = [int(row[1])]
                region_ends = [int(row[2])]
            aln = mafs_indexed[chrom].get_spliced(region_starts, region_ends, strand)
            aln_sp = {}
            for i in range(len(aln._records)):
                sp = aln._records[i].name.split('.')[0]
                aln._records[i].name = '{}.{}'.format(sp, seq_name)
                aln._records[i].id = '{}.{}'.format(sp, seq_name)
                aln_sp[sp] = i
            rec_new = list()
            if splist != []:
                for i in splist:
                    if i in aln_sp:
                        istr = FastaIO.as_fasta_2line(aln._records[aln_sp[i]])
                        print(istr, end='', file=outh)
            else:  # when splist is empty, output the records as is
                for i in aln._records:
                    print(FastaIO.as_fasta_2line(i), end='', file=outh)
            print('', file=outh)
    outh.close()
    return 0


if __name__ == '__main__':
    mafs_indexed = load_maf_index(maf_path=sys.argv[1], refsp=sys.argv[2])
    splist = load_splist(maf_path=sys.argv[1])
    parse_bed(mafs_indexed, splist, *sys.argv[3:])

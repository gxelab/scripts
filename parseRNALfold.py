# -*- coding: utf-8 -*-
"""
parse the output of RNALfold

output:
seq_id local_struct_start local_struct_len local_struct_energy z-score
-------------------------------------------------------------------------------
@author: zh (mt1022)
@date: Sun Dec  6 19:22:15 2015
"""
import sys
import re


def RNALfold_iter(lfold_output):
    """
    """
    struct_line = re.compile(r'^(.*?)\s+\(\s?(\S+)\)\s+(\d+)\s+z=\s+(\S+)')
    rna_id = ''
    local_struct = list()
    for line in open(lfold_output, 'r'):
        if line.startswith('>'):
            rna_id = line.strip().split(' ')[0][1:]
        elif line.startswith('.'):
            m = struct_line.search(line.rstrip())
            struct = m.group(1)
            kcal_mol = m.group(2)
            start = m.group(3)
            z_score = m.group(4)
            s = [start, str(len(struct)), kcal_mol, z_score, struct]
            local_struct.append(s)
        elif line.startswith(' '):
            yield rna_id, local_struct
            local_struct = list()


def parse_RNALfold(in_file):
    """
    """
    for rna_id, local_struct in RNALfold_iter(in_file):
        for s in local_struct:
            print('\t'.join([rna_id] + s))
        pass
    return


if __name__ == '__main__':
    parse_RNALfold(*sys.argv[1:])

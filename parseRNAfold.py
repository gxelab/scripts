# -*- coding: utf-8 -*-
"""
parse RNAfold output
------------------------------------------------
@author: zh (mt1022)
@date: Sat Dec 12 14:46:46 2015
"""
import sys
import re
import fileinput


def rnafold_iter(file_list):
    """
    """
    seq_id = ''
    seq_struct = ''
    seq_mfe = ''
    struct = re.compile('^([\.+\(\)]+)\s+\(\s{0,}(\S*?)\)$')
    for line in fileinput.input(file_list):
        line = line.rstrip()
        if line.startswith('>'):
            if seq_id != '':
                yield seq_id, seq_mfe, seq_struct
                seq_id = ''
                seq_struct = ''
                seq_mfe = ''
            seq_id = line.split(' ')[0][1:]
        elif line[0] in '.+()':
            m = struct.search(line)
            seq_struct = m.group(1)
            seq_mfe = m.group(2)
    yield seq_id, seq_mfe, seq_struct


if __name__ == '__main__':
    g4 = re.compile('\+\+\+')
    hp = re.compile('^[\.\+\(]+\.[\.\+\)]+$')
    bare = re.compile('^\.+$')
    for record in rnafold_iter(sys.argv[1:]):
        if g4.search(record[2]):
            is_g4 = '{}_{}'.format(record[2].find('+') + 1, record[2].rfind('+') + 1)
        else:
            is_g4 = 'NA'
        # is_g4 = 'g4' if g4.search(record[2]) else 'NA'
        if hp.search(record[2]) and not bare.search(record[2]):
            is_hairpin = 'hp'
        else:
            is_hairpin = 'NA'
        res = list(record) + [is_g4, is_hairpin]
        print('\t'.join(res))

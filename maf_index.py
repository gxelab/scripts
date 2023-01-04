#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from Bio.AlignIO import MafIO


def maf_index(chr_maf, target_name):
    """
    index a single maf file.
    
    example 1: index a single file
        python maf_index.py chr2L.maf dm6.chr2L
    example 2: index all mafs in current dir
        for i in *.maf; do echo python maf_index.py $i dm6.${i%%.maf}; done >00_build_index.sh
    """
    maf_index = chr_maf + 'idx'
    idx = MafIO.MafIndex(maf_index, chr_maf, target_name)
    return 0


if __name__ == '__main__':
    maf_index(*sys.argv[1:])
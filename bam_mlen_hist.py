#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute the histogram of mapped length of reads from input bam

To only count reads from certain genomic intervals, please use the "-L"
parameter of samtools to supply a bed file. 
------------------------------------------------------------------------------
Created on Mon Nov 16 20:12:51 2020
@author: zh (mt1022)
"""
import sys
import pysam
from collections import Counter


def main(bam_path):
    """
    generate histogram of mapped read lengths
    """
    hist = Counter()
    bam = pysam.AlignmentFile(bam_path, "rb")
    
    for read in bam:
        if read.mapping_quality >= 10 and read.is_supplementary == False:
            hist[read.query_alignment_length] += 1
        
    for i in range(min(hist.keys()), max(hist.keys()) + 1):
        print("{}\t{}".format(i, hist[i]))
    return


if __name__ == '__main__':
    main(*sys.argv[1:])

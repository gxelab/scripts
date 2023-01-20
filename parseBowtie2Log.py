# -*- coding: utf-8 -*-
"""
parse bowtie log file
------------------------------------------------
@author: zh (mt1022)
@date: Tue Aug 25 21:27:53 2015
"""
import sys
import re

bowtie_log = open(sys.argv[1], 'r')
infile = re.sub('.*/', '', sys.argv[1])

# extract total reads
total_reads = next(bowtie_log).rstrip().split(' ')[0]

# extract unmapped reads and percent
next(bowtie_log)
un_reads, un_percent = next(bowtie_log).rstrip().split()[:2]
un_percent = un_percent.translate(str.maketrans('', '', '()%'))

# print results
res = [infile, total_reads, un_reads, un_percent]
print('\t'.join(res))

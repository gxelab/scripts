# -*- coding: utf-8 -*-
"""
parse STAR log file
------------------------------------------------
@author: zh (mt1022)
@date: Wed Aug 26 10:37:33 2015
"""

import sys
import re

log_file = re.sub('.*/', '', sys.argv[1])
for line in open(sys.argv[1], 'r'):
    m = re.search('Number of input reads \|\t(\d+)', line)
    if m:
        input_reads = m.group(1)
    m = re.search('Uniquely mapped reads number \|\t(\d+)', line)
    if m:
        uniq_map = m.group(1)
    m = re.search('Uniquely mapped reads % \|\t([\d\.]+)', line)
    if m:
        uniq_percent = m.group(1)
    m = re.search('Number of reads mapped to multiple loci \|\t(\d+)', line)
    if m:
        multi_map = m.group(1)
    m = re.search('% of reads mapped to multiple loci \|\t([\d\.]+)', line)
    if m:
        multi_percent = m.group(1)
        break

s = [log_file, input_reads, uniq_map, uniq_percent, multi_map, multi_percent]
print('\t'.join(s))

# -*- coding: utf-8 -*-
"""
parse cutadapt log file
------------------------------------------------
@author: zh (mt1022)
@date: Tue Aug 25 21:52:47 2015
"""
import sys
import re

log_file = re.sub('.*/', '', sys.argv[1])
for line in open(sys.argv[1], 'r'):
    m = re.search('Total reads processed: *([\d,]+)', line)
    if m:
        total_reads = m.group(1).replace(',', '')
    m = re.search('Reads with adapters: *([\d,]+) \((\d+\.\d+)%\)', line)
    if m:
        with_adapt = m.group(1).replace(',', '')
        with_adapt_percent = m.group(2)
    m = re.search('Reads written \(passing filters\): *([\d,]+) \((\d+\.\d+)%\)', line)
    if m:
        reads_wirtten = m.group(1).replace(',', '')
        reads_written_percent = m.group(2)
        break

s = [log_file, total_reads, with_adapt, with_adapt_percent,
     reads_wirtten, reads_written_percent]
print("\t".join(s))

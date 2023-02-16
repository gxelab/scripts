import sys
import re
import pyBigWig
import pandas as pd
import numpy as np


def main(bw_forward_path, bw_reverse_path, bed_path):
    """
    Note:
        pd will read chromosomes names into numbers depending on the values
        in chroms in the columns, while bw.chroms().keys() is str. This has
        been updated (notified by WuCC)
    """
    bw_fw = pyBigWig.open(bw_forward_path)
    bw_rc = pyBigWig.open(bw_reverse_path)
    bed = pd.read_table(bed_path, header=None)

    for index, row in bed.iterrows():
        block_size = np.array([int(i) for i in row[10].split(',')[:-1]])
        block_start = np.array([int(i) for i in row[11].split(',')[:-1]])
        intv_start = row[1] + block_start
        intv_end = intv_start + block_size

        if row[5] == '+':
            bw = bw_fw
        else:
            bw = bw_rc
        chroms = bw.chroms().keys()

        chrom = str(row[0])  # chrom = re.sub(r'^chr', '', row[0])
        if chrom in chroms:
            ntcov = np.empty((0,))
            for i, j in zip(intv_start, intv_end):
                values = np.array(bw.values(chrom, i, j))
                values[np.isnan(values)] = 0
                ntcov = np.concatenate((ntcov, values))
            if row[5] == '-':
                ntcov = np.flip(ntcov)
        else:
            ntcov = np.zeros(sum(block_size))

        ntcov_char = ','.join(np.char.mod('%f', ntcov))
        print('\t'.join(str(i) for i in row[:6]), ntcov.size, ntcov_char, sep='\t')

if __name__ == '__main__':
    main(*sys.argv[1:])

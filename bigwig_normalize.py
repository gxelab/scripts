import sys
import numpy as np
import pyBigWig


def bw_write_chrom(bw, chrom, cov):
    """write the coverage of a chromosome into the bw file"""
    # note: np.array supports slice/modify by index vector
    run_end = np.append(np.where(np.diff(cov) != 0), len(cov)-1) + 1  # 1-based ends
    run_start = run_end - np.diff(np.append(0, run_end))  # 0-based starts
    run_value = cov[run_start]
    # ignore 0 values
    non_zero  = run_value > 0
    run_start = run_start[non_zero]
    run_end   = run_end[non_zero]
    run_value = run_value[non_zero]

    if len(run_value) > 0:
        bw.addEntries([chrom] * np.sum(non_zero),
                      run_start, ends = run_end, values = run_value)
    return


def main(bw_in_path, bw_out_path):
    bw_in = pyBigWig.open(bw_in_path)
    bw_out = pyBigWig.open(bw_out_path, 'w')
    # write header
    bw_header = [(i, bw_in.chroms(i)) for i in bw_in.chroms().keys()]
    bw_out.addHeader(bw_header)

    # count total reads
    total = 0
    for i, j in bw_header:
        total += np.nansum(bw_in.values(i, 0, j))

    # write normalized values
    for i, j in bw_header:
        new_values = bw_in.values(i, 0, j)
        new_values = np.nan_to_num(new_values)/total*1000000
        bw_write_chrom(bw_out, i, new_values)
    bw_out.close()
    return


if __name__ == '__main__':
    main(*sys.argv[1:])

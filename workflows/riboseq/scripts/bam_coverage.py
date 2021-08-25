"""

----------------------------------------
@author: zh (mt1022)
@date: Fri Nov 01 2019
"""

import sys
import argparse
import numpy as np
import pysam
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


def main(bam_path, output_prefix, fixed_offset,
             psite_offset_file, rlen_lower=27, rlen_upper=34):
    """
    calculate the coverage for plus strand and minus strand seperately
    """
    # try to open bam file
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    # create bw file and add header

    bw_fw = pyBigWig.open('{}_fw.bw'.format(output_prefix), 'w')
    bw_rc = pyBigWig.open('{}_rc.bw'.format(output_prefix), 'w')
    bw_header = [(i, bamfile.get_reference_length(i)) for i in bamfile.references]
    bw_fw.addHeader(bw_header)
    bw_rc.addHeader(bw_header)

    # parse P-site offset file
    psite_offset = dict()
    if psite_offset_file is None:
        for i in range(rlen_lower, rlen_upper + 1):
            psite_offset[i] = fixed_offset
    else:
        with open(psite_offset_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                ary = line.rstrip().split()
                try:
                    psite_offset[int(ary[0])] = int(ary[1])
                except ValueError:
                    continue

    # processing the bam file
    refname_curr = None  # name of current reference
    for read in bamfile:
        # mapping quality and read length QC
        if read.mapping_quality < 10 or read.is_supplementary:
            continue
        if read.query_alignment_length < rlen_lower or read.query_alignment_length > rlen_upper:
            continue

        # write bw before processing a new chromosome
        if read.reference_name != refname_curr:
            if refname_curr != None:
                bw_write_chrom(bw_fw, refname_curr, cov_plus)
                bw_write_chrom(bw_rc, refname_curr, cov_minus)
            refname_curr = read.reference_name
            cov_plus = np.zeros(bamfile.get_reference_length(read.reference_name))
            cov_minus = np.zeros(bamfile.get_reference_length(read.reference_name))

        # get P-site offset for read of this length
        # offset = 1 means count coverage for the 2nd nucleotide, ignore the first one
        # offset = 0 means no offset and count the 5' most position (helpful for writing next few lines)
        offset = psite_offset[read.query_alignment_length]
        # store coverage of current read to the coverage vectors
        # coverage vectors are 0-based
        # get_reference_positions() return 0-based positions on the reference
        # (https://pysam.readthedocs.io/en/latest/api.html#introduction)
        if read.is_reverse:
            cov_minus[ read.get_reference_positions()[-(offset + 1)] ] += 1
        else:
            cov_plus[ read.get_reference_positions()[offset] ] += 1

    # write last chromosome if not empty (note: at least one of them is not empty)
    if np.sum(cov_plus) > 0:
        bw_write_chrom(bw_fw, read.reference_name, cov_plus)
    if np.sum(cov_minus) > 0:
        bw_write_chrom(bw_rc, read.reference_name, cov_minus)

    # close bw files
    bw_fw.close()
    bw_rc.close()

    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog = 'python bam_coverage.py',
        description = 'Calculate coverage for Ribo-Seq data',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # required arguments (input and output)
    parser.add_argument('bam',
                        metavar='bam',
                        help='path to input bam [bam header required]')
    parser.add_argument('output_prefix',
                        metavar='output_prefix',
                       help='prefix of output files [output_prefix_(fw|rc).bw]')
    # define P-site offset
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-p', '--poffset',
                       type=int,
                       help='fixed P-site offset, for example, 12')
    group.add_argument('-v', '--variable',
                       help='for variable p-site offset, path to the offset file ('
                       'The file should be tab/space-seperated with two columns of '
                       'read length and offset)')
    # optional arguments
    parser.add_argument('-l', '--lower',
                        type=int,
                        default=27,
                        help='lower bound for mapped read length')
    parser.add_argument('-u', '--upper',
                        type=int,
                        default=34,
                        help='upper bound for mapped read length')
    
    args = parser.parse_args()
    if args.lower > args.upper:
        sys.exit('Lower bound of read length is larger than upper bound!')
    main(
        bam_path         =args.bam,
        output_prefix    =args.output_prefix,
        fixed_offset     =args.poffset,
        psite_offset_file=args.variable,
        rlen_lower       =args.lower,
        rlen_upper       =args.upper)

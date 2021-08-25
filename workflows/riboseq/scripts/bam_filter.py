"""
Filter BAM/SAM files
----------------------------------------
@author: zh (mt1022)
@date: Wed Aug 25 2021
"""

import sys
import argparse
import pysam


def main(bam_in, sam_out = '-', bout = False, qual=10, rlen_lower=27, rlen_upper=34):
    """
    filter bam
    """
    if bam_in[-3:] == 'bam':
        bam = pysam.AlignmentFile(bam_in, "rb")
    else:
        bam = pysam.AlignmentFile(bam_in, 'r')
    if bout:
        out = pysam.AlignmentFile(sam_out, 'wb', header=bam.header)
    else:
        out = pysam.AlignmentFile(sam_out, 'w', header=bam.header)

    for align in bam:
        # mapping quality and read length QC
        if (not align.is_supplementary and align.mapping_quality >= qual and
                align.query_alignment_length >= rlen_lower and
                align.query_alignment_length <= rlen_upper):
            out.write(align) # print(align.tostring(bam))
    bam.close()
    out.close()
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog = 'python bam_filter.py',
        description = 'Filter BAM/SAM file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # required arguments (input and output)
    parser.add_argument('bam',
                        metavar='bam',
                        help='path to input bam')
    # optional arguments
    parser.add_argument('-o', '--output',
                        type=str,
                        default='-',
                        help='output file name')
    parser.add_argument('-b', '--bamout',
                        action='store_true',
                        help='output in bam format')
    parser.add_argument('-q', '--qual',
                        type=int,
                        default=0,
                        help='min mapping quality')
    parser.add_argument('-l', '--lower',
                        type=int,
                        default=1,
                        help='lower bound for mapped read length')
    parser.add_argument('-u', '--upper',
                        type=int,
                        default=1000,
                        help='upper bound for mapped read length')
    
    args = parser.parse_args()
    if args.lower > args.upper:
        sys.exit('Lower bound of read length is larger than upper bound!')
    main(
        bam_in      = args.bam,
        sam_out     = args.output,
        bout        = args.bamout,
        qual        = args.qual,
        rlen_lower  = args.lower,
        rlen_upper  = args.upper)

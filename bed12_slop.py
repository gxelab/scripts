"""
adjust the size of interval stored in bed12 format;
=================================
author: mt1022 (zh)
date: 2023-12-13 20:50
"""
import sys
import argparse


def slop(bed_stream, left=0, right=0, out_stream=sys.stdout):
    """equivalent to bedtools slop for bed12"""
    for line in bed_stream:
        ary = line.strip().split('\t')
        blocks = ary[9]
        blocks_size = [int(i) for i in ary[10].rstrip(',').split(',')]
        blocks_start = [int(i) for i in ary[11].rstrip(',').split(',')]
        if sum(blocks_size) <= right + left:
            continue  # discard too short reads
        region_start = int(ary[1])
        # make chrosomal intervals
        intervals = []
        for i, j in zip(blocks_size, blocks_start):
            intervals.append((region_start + j, region_start + j + i - 1))
        # pruning left end
        if left >= 0:
            l = left
            while l > 0:
                i, j = intervals[0]
                if j - i + 1 > l:
                    intervals[0] = (i + l, j)
                    l = 0
                else:
                    intervals = intervals[1:]
                    l -= (j - i + 1)
        else:
            i, j = intervals[0]
            intervals[0] = (i + left, j)
        # pruning right end
        if right >= 0:
            r = right
            while r > 0:
                i, j = intervals[-1]
                if j - i + 1 > r:
                    intervals[-1] = (i, j - r)
                    r = 0
                else:
                    intervals = intervals[:-1]
                    r -= (j - i + 1)
        else:
            i, j = intervals[-1]
            intervals[-1] = (i, j - right)
        # update intervals
        new_blocks_size = [j - i + 1 for i, j in intervals]
        new_region_start = intervals[0][0]
        new_region_end = intervals[-1][1] + 1
        new_blocks_start = [i - new_region_start for i, j in intervals]
        # print output
        out = [ary[0], new_region_start, new_region_end] + ary[3:6] + [
            new_region_start, new_region_end, '255,0,0'] + [
            len(intervals), ','.join(str(i) for i in new_blocks_size),
            ','.join(str(i) for i in new_blocks_start)
        ]
        out = [str(i) for i in out]
        print('\t'.join(out), file=out_stream)
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='bed12tools slop')
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('-l', '--left', type=int, default=0,
                        help='The number of base pairs to remove from left side')
    parser.add_argument('-r', '--right', type=int, default=0,
                        help='The number of base pairs to remove from right side')
    parser.add_argument('-o', '--output', nargs='?', type=argparse.FileType('w'), default=sys.stdout)

    args = parser.parse_args()
    slop(args.input, args.left, args.right, args.output)

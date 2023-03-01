# -*- coding: utf-8 -*-
# Author: LeiTY
import random
import pysam
from collections import Counter
import click


@click.command(context_settings=dict(help_option_names=['-h', '--help'], show_default=True))
@click.argument('bam_path', type=click.STRING)
@click.option('-p', '--prefix', default='split_tmp', type=click.STRING,
              help='output prefix')
def bam_split(bam_path, prefix = 'split_tmp'):
    """
    Random split a bam file into two bam files based on read names.
    
    The two output files have the same number of reads.
    All the alignment of a single read will be put into the same file.
    """
    
    file = pysam.AlignmentFile(bam_path, 'rb')

    out1 = f'{prefix}.1.bam'
    out2 = f'{prefix}.2.bam'

    outfile1 = pysam.AlignmentFile(out1, 'wb', template = file)
    outfile2 = pysam.AlignmentFile(out2, 'wb', template = file)

    name = Counter()

    for line in file:
        name[line.query_name] = 1

    select_name = random.sample(list(name.keys()), len(name)//2)

    for i in select_name:
        name[i] = 2


    file.reset()
    for line in file:
        if name[line.query_name] == 1:
            outfile1.write(line)
        else:
            outfile2.write(line)


    file.close()       
    outfile1.close()
    outfile2.close()
        
    return
    
    
if __name__ == '__main__':
    bam_split()
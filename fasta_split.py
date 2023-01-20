import sys, re, argparse


def read_fasta(iterable):
    block = []
    for line in iterable:
        if line[0] == '>':
            if block:
                yield block
            block = [line]
        else:
            block.append(line)
    if block:
        yield block


def split_fasta(fasta, seqs, prefix, digits):
    """
    Split fasta files into files with fixed number of sequences
    """
    if prefix is None:
        prefix = xpath
    ofile_formatter = prefix + '{:0' + str(digits) + 'd}' + '.fa'

    n_seq = 0
    n_file = 1
    for record in read_fasta(fasta):
        if n_seq == 0:
            ofile_now = open(ofile_formatter.format(n_file), 'wt')
        out = ''.join(record)
        print(out, file = ofile_now, end='')
        n_seq += 1
        if n_seq == seqs:
            ofile_now.close()
            n_file += 1
            n_seq = 0
    ofile_now.close()
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog = 'python fasta_split.py',
        description = 'split fasta to files with n seqs each',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('fasta', nargs='?',
                        type=argparse.FileType('r'), default=sys.stdin,
                        help='input fasta file')
    parser.add_argument('-n', '--number',
                        type=int, default=1,
                        help='number of sequences in each split')
    parser.add_argument('-p', '--prefix',
                        type=str, default='fasplit_',
                        help='output prefix (as in prefix001.fa, etc)')
    parser.add_argument('-d', '--digits',
                        type=int, default=3,
                        help='number of digits to append to each file')
    
    args = parser.parse_args()
    split_fasta(fasta=args.fasta, seqs=args.number,
                prefix=args.prefix, digits=args.digits)

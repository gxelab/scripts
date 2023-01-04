import sys
import fileinput
from re import split


def mfa_iter(mfa_file):
    """
    iterate over MAF derived multi-block fasta file
    """
    header = ''
    seqs = dict()
    with fileinput.input(mfa_file) as f:
        for line in f:
            if line.startswith('\n'):
                if header != '':
                    yield(header, seqs)
                header = ''
                seqs = dict()
            elif line.startswith('>'):
                tmp = split(' |\.', line.strip()[1:])
                src = tmp[0]
                if header == '':
                    header = tmp[1]
                seq = next(f).rstrip().upper()
                seqs[src] = seq.translate({ord('#'): ord('-')})
        if header != '':
            yield(header, seqs)



def split_mfa(mfa_file, out_path):
    """
    perform split
    """
    if not out_path.endswith('/'):
        out_path = out_path + '/'
    for header, block in mfa_iter(mfa_file):
        new_file = out_path + header + '.fa'
        fh = open(new_file, 'w')

        for sp in block:
            seq_wgap = block[sp]
            if len(seq_wgap.replace('-', '')) > 0:
                print(">" + sp, file=fh)
                print(seq_wgap, file=fh)
        fh.close()
    return

if __name__ == '__main__':
    split_mfa(*sys.argv[1:])

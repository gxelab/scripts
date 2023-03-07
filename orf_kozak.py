import re
import gzip
import click
from Bio import SeqIO
import pandas as pd

def smart_open(path):
    """open plain text or gzipped file"""
    if path[-2:] == 'gz':
        return gzip.open(path, 'rt')
    else:
        return click.open_file(path, 'rt')


def strip_version(tx_name):
    """Strip transcript version"""
    return re.sub('\.\d+$', '', tx_name)


def read_fasta(path, ignore_version=False):
    """Construct a dict of input fasta"""
    fa = dict()
    with smart_open(path) as f:
        for record in SeqIO.parse(f, 'fasta'):
            if ignore_version:
                record.id = strip_version(record.id)
            fa[record.id] = str(record.seq)
    return fa


@click.command(context_settings=dict(help_option_names=['-h', '--help'], show_default=True))
@click.argument('orf_table', type=click.STRING)
@click.option('-f', '--fas_path', type=click.STRING, default=['-'], multiple=True,
              help=('transcript sequences in fasta. can be specified multiple times if ' +
              'sequences are stored in separate files (e.g., cdna.fa and ncrna.fa)'))
@click.option('-i', '--ignore_txversion', is_flag=True, default=False,
              help='ignore transcript version in ".\d+" format')
def get_kozak(orf_table, fas_path, ignore_txversion=True):
    """
    Extract Kozak sequence for each ORF
    
    \b
    ORF_TABLE: table of ORFs (must have columns orf_id, tx_name, tstart)
    output: stdout
    """
    fas = dict()
    for i in fas_path:
        fas = fas | read_fasta(i, ignore_version=ignore_txversion)
    orfs = pd.read_table(orf_table)
    print('\t'.join(['orf_id', 'tx_name', 'kozak_seq', 'kozak_score']))
    for row in orfs.itertuples():
        try:
            kozak = ('-' * max(0, 7 - row.tstart) + 
                fas[row.tx_name][max(0, row.tstart - 7):(row.tstart + 3)])
            score = (kozak[3] in 'AG') + (kozak[9] == 'G')
        except KeyError:
            kozak = '##########'
            score = 'NA'
        print('\t'.join([row.orf_id, row.tx_name, kozak, str(score)]))
    return


if __name__ == '__main__':
    get_kozak()

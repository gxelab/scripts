import sys
import click
import pandas as pd


def drop_index_col(df):
    if 'Unnamed: 0' in df.columns:
        df = df.drop(columns=['Unnamed: 0'])
    else:
        df = df.copy()
    return df


@click.command(context_settings=dict(help_option_names=['-h', '--help'], show_default=True))
@click.argument('orfs_table', type=click.STRING)
@click.option('-q', '--quant_table', type=click.STRING, default=None,
              help='results of orf_quant.py')
@click.option('-f', '--floss_table', type=click.STRING, default=None,
              help='results of orf_floss_cutoff.R')
@click.option('-k', '--kozak_table', type=click.STRING, default=None,
              help='results of orf_kozak.py')
@click.option('-o', '--output', type=click.File('w'), default=sys.stdout,
              help='output path')
def merge_tables(orfs_table, quant_table, floss_table, kozak_table, output):
    """
    merge ORF tables
    
    \b
    ORFS_TABLE: result of orf_type.py
    """
    orfs = pd.read_table(orfs_table, dtype={'chrom': 'string'})
    orfs = orfs[['orf_id', 'tx_name']]
    if quant_table is not None:
        quant = pd.read_table(quant_table, dtype={'chrom': 'string'})
        quant = drop_index_col(quant)
        orfs = pd.merge(orfs, quant, on=['orf_id', 'tx_name'], how='left')
    if floss_table is not None:
        floss = pd.read_table(floss_table)
        orfs = pd.merge(orfs, floss, on=['orf_id', 'tx_name'], how='left')
    if kozak_table is not None:
        kozak = pd.read_table(kozak_table)
        orfs = pd.merge(orfs, kozak, on=['orf_id', 'tx_name'], how='left')
    orfs.to_csv(output, sep='\t', float_format='%g', index=False)
    return


if __name__ == '__main__':
    merge_tables()

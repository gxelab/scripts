import pyBigWig as pbw
import pandas as pd
import numpy as np
import click


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.argument('bw_forward_path', type=click.STRING)
@click.argument('bw_reverse_path', type=click.STRING)
@click.argument('bed_path', type=click.STRING)
def main(bw_forward_path, bw_reverse_path, bed_path):
    """
    Coverage stats for ORFs in the bed file.
    """
    bw_fw = pbw.open(bw_forward_path)
    bw_rc = pbw.open(bw_reverse_path)
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

        chrom = str(row[0])
        if chrom in chroms:
            ntcov = np.empty((0,))
            for i, j in zip(intv_start, intv_end):
                values = np.nan_to_num(bw.values(chrom, i, j))
                ntcov = np.concatenate((ntcov, values))
            if row[5] == '-':
                ntcov = np.flip(ntcov)
        else:
            ntcov = np.zeros(sum(block_size))

        # values after excluding the first and last 3 values (remove the first and last codon)
        ntcov_wos = ntcov[3:-3]

        print('\t'.join(str(i) for i in row[:6]), ntcov.size,
              # sum_all, sum_f0, sum_f1, sum_f2,
              ntcov.sum(), ntcov[::3].sum(), ntcov[1::3].sum(), ntcov[2::3].sum(),
              # sum_wos, sum_wos_f0, sum_wos_f1, sum_wos_f2
              ntcov_wos.sum(), ntcov_wos[::3].sum(), ntcov_wos[1::3].sum(), ntcov_wos[2::3].sum(), sep='\t')


if __name__ == '__main__':
    main()

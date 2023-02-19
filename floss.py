from collections import Counter
from multiprocessing import Pool
import click
import pysam
import numpy as np
import pandas as pd
import gppy.gtf as gtf



@click.command(context_settings=dict(help_option_names=['-h', '--help'], show_default=True))
@click.argument('path_bam', type=click.STRING)
@click.argument('path_gtf', type=click.STRING)
@click.argument('txinfo_table', type=click.STRING)
@click.argument('cds_bed', type=click.STRING)
@click.argument('orf_table', type=click.STRING)
@click.option('-o', '--output', type=click.STRING, default=None,
              help='output table')
@click.option('-p', '--n_cpus', type=click.INT, default=1,
              help='whether to exclude annotated CDSs and CDS isoforms')
def compute_floss(path_bam, path_gtf, txinfo_table, cds_bed, orf_table, output=None, n_cpus=1):
    """
    Compute Fragment Length Organization Similarity Score (FLOSS)

    \b
    path_bam: genomic aligment
    path_gtf: genome annotation
    txinfo_table: txinfo abtained with gppy txinfo -f -g
    cds_bed: cds genomic intervals merged in a strand-specific manner
    orf_table: ORF table generated by ncorf_classifier[version].py
    output: output results
    \b
    example command to generated required input tables:
    ```
    gppy convert2bed -t cds -g Drosophila_melanogaster.BDGP6.32.52.gtf \\
        | bedtools bed12tobed6 | bedtools sort | bedtools merge -s -c 4,5,6 -o first > cdsmerge.bed
    gppy txinfo -g Drosophila_melanogaster.BDGP6.32.52.gtf -f >Drosophila_melanogaster.BDGP6.32.52.gtf.txinfo
    ```

    Reference: Ingolia et al., 2014, Cell Reports
    """
    gannot = gtf.parse_gtf(path_gtf)

    orfs = pd.read_table(orf_table)
    orfs['givs'] = orfs.apply(lambda r: gannot[r.tx_name].tiv_to_giv(r.tstart, r.tend), axis=1)

    txinfo = pd.read_table(txinfo_table)
    txinfo_pcg = txinfo[(txinfo.transcript_biotype == 'protein_coding') & (txinfo.cds_len > 30)]
    txinfo_pcg = txinfo_pcg.assign(cds_start=txinfo.utr5_len + 1, cds_end=txinfo.utr5_len + txinfo.cds_len)
    txinfo_pcg['givs'] = txinfo_pcg.apply(lambda r: gannot[r.tx_name].tiv_to_giv(r.cds_start, r.cds_end), axis=1)
    txinfo_pcg = txinfo_pcg[['gene_id', 'tx_name', 'chrom', 'strand', 'givs']]

    bam = pysam.AlignmentFile(path_bam, 'rb')
    cds = pd.read_table(cds_bed, header=None)
    cds = cds.astype({0: 'string'})
    
    cds_dist = Counter()

    # reference read length frequency distribution
    def worker(row):
        rlen_dist = Counter()
        bam = pysam.AlignmentFile(path_bam, 'rb')
        for align in bam.fetch(contig=row[1], start=row[2], stop=row[3]):
            if align.is_reverse == (row[6] == '-'):
                rlen_dist[align.query_alignment_length] += 1
        return rlen_dist

    if n_cpus == 1:
        for row in cds.itertuples():
            for align in bam.fetch(contig=row[1], start=row[2], stop=row[3]):
                if align.is_reverse == (row[6] == '-'):
                    cds_dist[align.query_alignment_length] += 1
    else:
        with Pool(processes=n_cpus) as pool:
            for counter in pool.map(worker, cds.itertuples(name=None)):
                cds_dist += counter
    
    cds_dist = pd.DataFrame.from_dict(cds_dist, orient='index').reset_index()
    cds_dist.columns = ['len', 'cnt']
    cds_dist = cds_dist.sort_values(by='len')
    cds_dist['freq'] = cds_dist.cnt / cds_dist.cnt.sum()
    ref_lens = cds_dist.len.to_numpy()
    ref_freq = cds_dist.freq.to_numpy()
    
    # compute FLOSS
    def floss(row):
        orf_dist = Counter()
        for region in row.givs:
            # Region start/end is 1-based
            for align in bam.fetch(contig=row.chrom, start=region.start - 1, end=region.end):
                if align.is_reverse != (row.strand == '-'):
                    continue
                offset = align.get_tag('PS')
                if align.is_reverse:
                    psite = align.get_reference_positions()[-(offset + 1)]
                else:
                    psite = align.get_reference_positions()[offset]
                # aligned positions is 0-based
                if region.start - 1 <= psite < region.end:
                    orf_dist[align.query_alignment_length] += 1
        orf_dist = np.array([orf_dist[i] for i in ref_lens])
        psite_tot = orf_dist.sum()
        orf_dist = orf_dist / (psite_tot + 0.001)
        score = np.abs(orf_dist - ref_freq).sum()/2
        return psite_tot, score


    def floss_worker(row):
        """
        row: 1-orf_id, 2-tx_name, 3-chrom, 4-strand, 5-givs
        """
        orf_dist = Counter()
        bam = pysam.AlignmentFile(path_bam, 'rb')
        for region in row[5]:
            # Region start/end is 1-based
            for align in bam.fetch(contig=row[3], start=region.start - 1, end=region.end):
                if align.is_reverse != (row[4] == '-'):
                    continue
                offset = align.get_tag('PS')
                if align.is_reverse:
                    psite = align.get_reference_positions()[-(offset + 1)]
                else:
                    psite = align.get_reference_positions()[offset]
                # aligned positions is 0-based
                if region.start - 1 <= psite < region.end:
                    orf_dist[align.query_alignment_length] += 1
        orf_dist = np.array([orf_dist[i] for i in ref_lens])
        psite_tot = orf_dist.sum()
        orf_dist = orf_dist / (psite_tot + 0.001)  # pseudocount to avoid divid by zero
        score = np.abs(orf_dist - ref_freq).sum()/2
        return psite_tot, score
    
    orfs_short = orfs[['orf_id', 'tx_name', 'chrom', 'strand', 'givs']].copy()
    if n_cpus == 1:
        orfs_short[['floss_cnt', 'floss']] = orfs_short.apply(floss, axis=1, result_type='expand')
        txinfo_pcg[['floss_cnt', 'floss']] = txinfo_pcg.apply(floss, axis=1, result_type='expand')
    else:
        with Pool(processes=n_cpus) as pool:
            orfs_short[['floss_cnt', 'floss']] = pool.map(floss_worker, orfs_short.itertuples(name=None))
            txinfo_pcg[['floss_cnt', 'floss']] = pool.map(floss_worker, txinfo_pcg.itertuples(name=None))
    
    txinfo_pcg.rename(columns={'gene_id': 'id'}, inplace=True)
    orfs_short.rename(columns={'orf_id': 'id'}, inplace=True)

    floss_tab = pd.concat([orfs_short.assign(type='ORF'), txinfo_pcg.assign(type='CDS')], axis=0)
    floss_tab = floss_tab[['id', 'tx_name', 'floss_cnt', 'floss', 'type']]
    floss_tab.floss_cnt = floss_tab.floss_cnt.fillna(0)
    floss_tab = floss_tab.astype({'floss_cnt': 'int32'})
    floss_tab
    if output == None:
        output = f'{orf_table}.floss'
    floss_tab.to_csv(output, sep='\t', float_format='%g')
    return


if __name__ == '__main__':
    compute_floss()
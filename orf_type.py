import sys
import re
import gzip
from typing import List
from uuid import uuid4
from pathlib import Path
from os import path as os_path, remove as os_remove

import pandas as pd
import numpy as np
import pybedtools
from statsmodels.stats.multitest import multipletests
import click

import gtf as GTF


###############################################################################
# functions ###################################################################
###############################################################################
def get_exons_bed(gtf_file, outpath):
    """
    convert gtf to bed6 of exons
    """
    if gtf_file.endswith('.gz'):
        f = gzip.open(gtf_file, 'rt')
    elif gtf_file == '-':
        f = sys.stdin
    else:
        f = open(gtf_file)

    regex_gid = re.compile(r'gene_id "(.*?)"')
    regex_tid = re.compile(r'transcript_id "(.*?)"')
    with open(outpath, 'w') as bedh:
        for line in f:
            if line[0] == '#':
                continue
            ary = line.strip().split('\t')
            m_tid = regex_tid.search(ary[8])
            m_gid = regex_gid.search(ary[8])
            if ary[2] != 'exon':
                continue
            bstart = str(int(ary[3]) - 1)
            out = [ary[0], bstart, ary[4], m_tid.group(1), m_gid.group(1), ary[6]]
            print('\t'.join(out), file=bedh)
    f.close()
    return


def get_tiv_gene(row):
    tx = gtf[row.tx_name]
    return (tx.gpos_to_tpos(row.gstart)[0], tx.gpos_to_tpos(row.gend)[0],
            tx.gene.gene_id, tx.gene.gene_name)


def get_giv_gene(row):
    tx = gtf[row.tx_name]
    return (tx.tpos_to_gpos(row.tstart), tx.tpos_to_gpos(row.tend),
            tx.gene.gene_id, tx.gene.gene_name)


def read_riborf(path):
    """Read RibORF output"""    
    # split original ID
    riborf = pd.read_table(path, dtype={'chrom': str})
    riborf.rename(columns={'orfID': 'orf_id', 'pred.pvalue': 'p'}, inplace=True)
    riborf[['tx_name', 'tiv', 'orf_type_ori', 'start_codon']] = riborf.orf_id.str.split('|', expand=True)[[0, 2, 3, 4]]
    riborf['tx_name'] = riborf.tx_name.str.split(':', expand=True)[0]

    # add giv, tiv, gene info
    riborf[['tstart', 'tend']] = riborf.tiv.str.split(':', expand=True)[[1, 2]].astype('int')
    riborf['tend'] = riborf.tend -1
    riborf[['gstart', 'gend', 'gene_id', 'gene_name']] = riborf.apply(get_giv_gene, axis=1, result_type='expand')
    
    riborf = riborf[['orf_id', 'tx_name', 'tstart', 'tend', 'chrom', 'gstart', 'gend', 'strand',
                     'start_codon', 'orf_type_ori', 'gene_id', 'gene_name', 'p']]
    q10 = np.quantile(riborf[riborf.orf_type_ori == 'canonical'].p, q = 0.1)
    riborf = riborf[riborf.p > q10]
    riborf.drop(columns=['p'], inplace=True)
    return riborf


def read_ribocode(path):
    """
    Read RIBOCODE results
    
    note: adjusted pval in collapsed output is meaningless since it's calculated after
    excluding entries with raw p value > 0.05.
    """
    ribocode = pd.read_table(path, type={'chrom': 'string'})
    if 'start_codon' not in ribocode.columns:
        ribocode['start_codon'] = 'ATG'
    ribocode = ribocode[['ORF_ID', 'transcript_id', 'ORF_tstart', 'ORF_tstop', 'chrom', 'ORF_gstart',
                         'ORF_gstop', 'strand', 'start_codon', 'ORF_type', 'gene_id', 'gene_name']]
    ribocode.rename(columns={
        'ORF_ID': 'orf_id', 'transcript_id': 'tx_name', 'ORF_tstart': 'tstart', 'ORF_tstop': 'tend',
        'ORF_gstart': 'gstart', 'ORF_gstop': 'gend', 'ORF_type': 'orf_type_ori'}, inplace=True)
    return ribocode


def get_giv_boundary(giv_str):
    starts, ends = zip(*[i.split('-') for i in giv_str.split('|')])
    starts = list(map(int, starts))
    ends = list(map(int, ends))
    return min(starts) + 1, max(ends)


def read_price(path):
    """Read PRICE results"""
    price = pd.read_table(path)
    price['tx_name'] = price.Id.str.split('_', expand=True)[0]
    price[['chrom_strand', 'givs']] = price.Location.str.split(':', expand=True)
    price['chrom'] = price.chrom_strand.str[:-1]
    price['strand'] = price.chrom_strand.str[-1]
    price[['giv_start', 'giv_end']] = price.givs.apply(get_giv_boundary).tolist()
    price['gstart'] = np.where(price.strand == '+', price.giv_start, price.giv_end)
    price['gend'] = np.where(price.strand == '+', price.giv_end, price.giv_start)
    price[['tstart', 'tend', 'gene_id', 'gene_name']] = price.apply(get_tiv_gene, axis=1, result_type='expand')
    price['padj'] = multipletests(price['p value'], method='fdr_bh')[1]
    price = price[price.padj < 0.05]
    price.rename(columns={'Id': 'orf_id', 'Codon': 'start_codon', 'Type': 'orf_type_ori'}, inplace=True)
    price = price[['orf_id', 'tx_name', 'tstart', 'tend', 'chrom', 'gstart', 'gend',
                   'strand', 'start_codon', 'orf_type_ori', 'gene_id', 'gene_name']]
    return price


def read_ribotish(path):
    """Read Ribo-TISH results"""
    ribotish = pd.read_table(path)
    ribotish[['chrom', 'giv', 'strand']] = ribotish.GenomePos.str.split(':', expand=True)
    ribotish[['giv_start', 'giv_end']] = ribotish.giv.str.split('-', expand=True).astype(int)
    ribotish['gstart'] = np.where(ribotish.strand == '+', ribotish.giv_start + 1, ribotish.giv_end)
    ribotish['gend'] = np.where(ribotish.strand == '+', ribotish.giv_end, ribotish.giv_start + 1)
    ribotish.rename(columns={'Tid': 'tx_name', 'StartCodon': 'start_codon', 'TisType': 'orf_type_ori'}, inplace=True)
    ribotish[['tstart', 'tend', 'gene_id', 'gene_name']] = ribotish.apply(get_tiv_gene, axis=1, result_type='expand')
    ribotish['orf_id'] = ribotish.tx_name + '_' + ribotish.tstart.astype(str) + '_' + ribotish.tend.astype(str)
    ribotish = ribotish[[
        'orf_id', 'tx_name', 'tstart', 'tend', 'chrom', 'gstart', 'gend',
        'strand', 'start_codon', 'orf_type_ori', 'gene_id', 'gene_name']]
    return ribotish


def read_ribotricer(path):
    """Read RiboTRicer results"""
    ribotricer = pd.read_table(path, dtype={'chrom': 'string'})
    ribotricer.drop(columns = 'profile', inplace=True)
    ribotricer.rename(columns={'ORF_ID': 'orf_id', 'transcript_id': 'tx_name', 'ORF_type': 'orf_type_ori'}, inplace=True)
    ribotricer[['giv_start', 'giv_end']] = ribotricer.orf_id.str.split('_', expand = True).iloc[:,1:3].astype(int)
    ribotricer['gstart'] = np.where(ribotricer.strand == '+', ribotricer.giv_start, ribotricer.giv_end)
    ribotricer['gend'] = np.where(ribotricer.strand == '+', ribotricer.giv_end, ribotricer.giv_start)
    ribotricer[['tstart', 'tend', 'gene_id', 'gene_name']] = ribotricer.apply(get_tiv_gene, axis=1, result_type='expand')
    ribotricer['tend'] = ribotricer.tend + 3
    ribotricer[['gstart', 'gend', 'gene_id', 'gene_name']] = ribotricer.apply(get_giv_gene, axis=1, result_type='expand')
    ribotricer = ribotricer[[
        'orf_id', 'tx_name', 'tstart', 'tend', 'chrom', 'gstart', 'gend',
        'strand', 'start_codon', 'orf_type_ori', 'gene_id', 'gene_name']]
    return ribotricer


def orf_typing(row):
    if row.txtype != 'mRNA':
        return f'{row.txtype}-ORF'
    else:
        if row.tstart < row.cds_start:
            if row.tend < row.cds_start:
                return 'uORF'
            elif row.tend < row.cds_end:
                return 'uoORF'
            elif row.tend == row.cds_end:
                return 'N_extension'
            else:  # row.tend > row.cds_end
                if row.phase_start == 0 or row.phase_end == 0:
                    return 'wCDS'
                else:
                    return 'wORF'
        elif row.tstart == row.cds_start:
            if row.tend < row.cds_end:
                return 'C_truncation'
            elif row.tend == row.cds_end:
                return 'CDS'
            else:  # row.tend > row.cds_end
                return 'C_extension'
        elif row.tstart <= row.cds_end:
            if row.tend < row.cds_end:
                if row.phase_start == 0 or row.phase_end == 0:
                    return 'iCDS'
                else:
                    return 'iORF'
            elif row.tend == row.cds_end:
                return 'N_truncation'
            else:  # row.tend > row.cds_end
                if row.phase_start == 0 or row.phase_end == 0:
                    return 'sCDS'
                else:
                    return 'doORF'
        else:  #row.tstart > row.cds_end
            return 'dORF'


# https://www.gencodegenes.org/pages/biotypes.html
# https://asia.ensembl.org/info/genome/genebuild/biotypes.html
BIOTYPES = dict(
     mRNA=['protein_coding'],
    ncRNA=[
        'Mt_rRNA', 'Mt_tRNA', 'miRNA', 'pre_miRNA', 'misc_RNA', 'rRNA', 'known_ncrna',
        'ribozyme', 'sRNA', 'scRNA', 'scaRNA', 'snRNA', 'snoRNA', 'ncRNA', 'tRNA', 'vault_RNA', 'vaultRNA'],
    lncRNA=[
        'lncRNA', '3prime_overlapping_ncRNA', 'antisense', 'antisense_RNA',
        'lincRNA', 'macro_lncRNA', 'processed_transcript', 'sense_intronic',
        'sense_overlapping', 'non_coding', 'bidirectional_promoter_lncRNA'],
    pseudogene=[
        'polymorphic_pseudogene', 'processed_pseudogene', 'pseudogene', 'retrotransposed',
        'transcribed_processed_pseudogene', 'transcribed_unitary_pseudogene',
        'transcribed_unprocessed_pseudogene', 'translated_processed_pseudogene',
        'translated_unprocessed_pseudogene', 'unitary_pseudogene', 'unprocessed_pseudogene',
        'Mt_tRNA_pseudogene', 'tRNA_pseudogene', 'snoRNA_pseudogene', 'snRNA_pseudogene',
        'scRNA_pseudogene', 'rRNA_pseudogene', 'misc_RNA_pseudogene', 'miRNA_pseudogene'],
    IGTR=[
        'IG_C_gene', 'IG_C_pseudogene', 'IG_D_gene', 'IG_J_gene', 'IG_J_pseudogene', 'IG_LV_gene',
        'IG_V_gene', 'IG_V_pseudogene', 'IG_pseudogene', 'TR_C_gene', 'TR_D_gene',
        'TR_J_gene', 'TR_J_pseudogene', 'TR_V_gene', 'TR_V_pseudogene'],
    misc=[
        'non_stop_decay', 'nonsense_mediated_decay', 'TEC', 'artifact',
        'disrupted_domain', 'ambiguous_orf', 'retained_intron', 'Readthrough',
        'protein_coding_LoF', 'protein_coding_CDS_not_defined'] 
    # biotypes not present in this list will be assigned to "misc" category
)


def get_biotype_table():
    biotype_table = dict(txtype=[], biotype=[])
    for k, v in BIOTYPES.items():
        biotype_table['txtype'].extend([k] * len(v))
        biotype_table['biotype'].extend(v)
    return pd.DataFrame(biotype_table)


CHROM_RM = ['MT', 'mitochondrion_genome']

###############################################################################
# main analysis ###############################################################
###############################################################################
@click.command(context_settings=dict(help_option_names=['-h', '--help'], show_default=True))
@click.argument('path_orf_pred', type=click.STRING)
@click.argument('path_gtf', type=click.STRING)
@click.argument('path_txinfo', type=click.STRING)
@click.option('-m', '--pred_method', default='bed12',
              type=click.Choice(['price', 'ribocode', 'riborf', 'ribotish', 'ribotricer', 'bed12', 'tabular']),
              help='ORF prediction method')
@click.option('-p', '--prefix', default='orfpred', type=click.STRING,
              help='output prefix')
@click.option('-c', '--exclude_cds', is_flag=True, default=False,
              help='exclude annotated CDSs and CDS isoforms')
@click.option('-u', '--unique_by_coord', is_flag=True, default=False,
              help='make input ORFs unique by genomic coordinates')
def main(path_orf_pred, path_gtf, path_txinfo, pred_method, prefix, exclude_cds=False, unique_by_coord=False):
    """
    A pipeline for the classification of non-canonical ORFs

    \b
    path_orf_pred: path to predicted ORFs
    path_gtf: path to gtf file
    path_txinfo: path to transcript info file (e.g., generated by extract_txinfo_ensembl.R)
    """
    # input assertion
    assert Path(path_gtf).is_file()
    assert Path(path_txinfo).is_file()
    assert pred_method in ['ribocode', 'ribotish', 'price', 'riborf', 'ribotricer', 'bed12', 'tabular']

    if pred_method == 'riborf':
        if Path(path_orf_pred).is_dir():
            path_orf_pred = os_path.join(path_orf_pred, 'repre.valid.pred.pvalue.parameters.txt')

    assert Path(path_orf_pred).is_file()

    # output
    path_output = prefix + '_processed.tsv'
    path_outraw = prefix + '_raw.tsv'

    # load annotation in GTF format ===============================================
    global gtf
    gtf = GTF.parse_gtf(path_gtf)

    # load transcript info ========================================================
    print('...load txinfo', file=sys.stderr)
    txinfo = pd.read_csv(path_txinfo, sep='\t')
    txinfo['tx_complete'] = ((txinfo.cds_start_nf | txinfo.cds_end_nf |
        txinfo.mrna_start_nf | txinfo.mrna_end_nf) != True)

    # extract ends of annotated CDS
    txinfo['cds_start'] = np.where(
        txinfo.transcript_biotype == 'protein_coding',
        txinfo.utr5_len + 1, np.nan)
    txinfo['cds_end'] = np.where(
        txinfo.transcript_biotype == 'protein_coding',
        txinfo.utr5_len + txinfo.cds_len, np.nan)

    cds_ends = pd.concat([
        (txinfo
        .query('transcript_biotype == "protein_coding" and not cds_start_nf')
        .assign(
            posn = txinfo.cds_start,
            type = 'cds_start')[['tx_name', 'posn', 'type', 'chrom']]),
        (txinfo
        .query('transcript_biotype == "protein_coding" and not cds_end_nf')
        .assign(
            posn = txinfo.cds_end,
            type = 'cds_end')[['tx_name', 'posn', 'type', 'chrom']]) ])
    cds_ends['gposn'] = cds_ends.apply(lambda r: gtf[r.tx_name].tpos_to_gpos(int(r.posn)), axis=1)

    # clean transcripts used to annotate ORF types
    # biotypes_keep = BIOTYPES['mRNA'] + BIOTYPES['ncRNA'] + BIOTYPES['lncRNA'] + BIOTYPES['pseudogene']
    # txinfo_clean = txinfo_raw[txinfo_raw.transcript_biotype.isin(biotypes_keep)]

    # read ORF prediction =========================================================
    print('...read ORFs', file=sys.stderr)
    if pred_method == 'ribocode':
        orfs_pred = read_ribocode(path_orf_pred)
    elif pred_method == 'riborf':
        orfs_pred = read_riborf(path_orf_pred)
    elif pred_method == 'price':
        orfs_pred = read_price(path_orf_pred)
    elif pred_method == 'ribotish':
        orfs_pred = read_ribotish(path_orf_pred)
    elif pred_method == 'ribotricer':
        orfs_pred = read_ribotricer(path_orf_pred)
    elif pred_method == 'tabular':
        orfs_pred = pd.read_csv(path_orf_pred, sep='\t')
    elif pred_method == 'bed12':
        orfs_bed12 = pybedtools.BedTool(path_orf_pred).filter(lambda item: item.chrom not in CHROM_RM).saveas()
        orfs_pred = orfs_bed12.to_dataframe(names=[
            'orf_chrom', 'orf_start', 'orf_end', 'orf_id', 'orf_score', 'orf_strand', 'thick_start',
            'thick_end', 'item_rgb', 'block_count', 'block_sizes', 'block_starts'])
        if orfs_pred.orf_id.duplicated().sum() > 0:
            sys.exit('The name (4th) column has duplicated values!')

    # save raw results
    orfs_pred.to_csv(path_outraw, sep='\t', index=False)

    # read ORF prediction =========================================================
    print('...convert to bed6', file=sys.stderr)
    if pred_method == 'bed12':
        # calc ORF length
        if orfs_pred['block_sizes'][0][-1] == ',':
            orfs_pred['orf_len'] = orfs_pred['block_sizes'].apply(lambda x: sum(map(int, x.split(',')[:-1])))
        else:
            orfs_pred['orf_len'] = orfs_pred['block_sizes'].apply(lambda x: sum(map(int, x.split(','))))
        # remove duplicated ORFs that have same genomic coordinates as previous ones
        if unique_by_coord:
            orfs_pred = orfs_pred.drop_duplicates(
                subset=['orf_chrom', 'orf_start', 'orf_end', 'orf_len'], ignore_index=True)
        # get ORF giv bed (6 columns)
        orfs_bed = orfs_bed12.bed12tobed6()
    else:
        # calc ORF length
        orfs_pred['orf_len'] = orfs_pred.tend - orfs_pred.tstart + 1
        # rmeove genes encoded by mitochrondria, which use non-standard genetic code
        orfs_pred = orfs_pred[~orfs_pred.chrom.isin(CHROM_RM)]
        # remove duplicated ORFs that have same genomic coordinates as previous ones
        if unique_by_coord:
            orfs_pred = orfs_pred.drop_duplicates(
                subset=['chrom', 'gstart', 'gend', 'orf_len'], ignore_index=True)

        # get ORF giv
        orfs_giv = orfs_pred.copy()
        orfs_giv['givs'] = orfs_giv.apply(lambda r: gtf[r.tx_name].tiv_to_giv(r.tstart, r.tend), axis=1)
        orfs_giv = orfs_giv.explode('givs', ignore_index=True)
        orfs_giv['giv_start'] = [r.start for r in orfs_giv.givs]
        orfs_giv['giv_end'] = [r.end for r in orfs_giv.givs]
        orfs_giv.drop(columns=['givs'], inplace=True)

        # convert ORF giv to a tempary bed file
        orfs_bed = orfs_giv.loc[:, ['chrom', 'giv_start', 'giv_end', 'orf_id', 'tx_name', 'strand']]
        orfs_bed["giv_start"] = orfs_bed.giv_start - 1
        orfs_bed = orfs_bed.astype({'giv_start': 'int64', 'giv_end': 'int64'})

        orfs_bed_path = '/tmp/' + str(uuid4()) + '.bed'
        orfs_bed.to_csv(orfs_bed_path, sep = '\t', header=False, index=False)
        orfs_bed = pybedtools.BedTool(orfs_bed_path)

    # remap ORFs ==================================================================
    print('...fetch compatible tx isoforms', file=sys.stderr)

    # convert transcript exons to a tempary bed file
    exons_bed = '/tmp/' + str(uuid4()) + '.bed'
    get_exons_bed(path_gtf, exons_bed)

    # pybedtools intersect
    orfs_intersect = orfs_bed.intersect(exons_bed, wao=True, s=True).to_dataframe(names=[
        'orf_chrom', 'orf_start', 'orf_end', 'orf_id', 'orf_score', 'orf_strand', 'exon_chrom',
        'exon_start', 'exon_end', 'exon_name', 'exon_score', 'exon_strand', 'overlap'],
        dtype={0: str, 6: str})

    # remove tempory files
    os_remove(exons_bed)
    if pred_method != 'bed12':
        os_remove(orfs_bed_path)

    # filtering overlap
    ## ORF intervals should be fully contained in exons
    orfs_intersect = orfs_intersect[orfs_intersect.orf_end - orfs_intersect.orf_start == orfs_intersect.overlap]

    ## aggregate overlaps
    intersect_fine = (orfs_intersect
        .groupby(by=['orf_id', 'exon_name'])
        .agg({'overlap': 'sum', 'orf_start': lambda x: x.min() + 1, 'orf_chrom': lambda x: x.iloc[0],
            'orf_strand': lambda x: x.iloc[0], 'orf_end': 'max'})
        .reset_index())
    intersect_fine = intersect_fine.rename(columns={
        'exon_name': 'tx_name', 'overlap': 'overlap_tot', 'orf_chrom': 'chrom', 'orf_strand': 'strand'})
    intersect_fine['overlap_start'] = np.where(intersect_fine.strand == '+', intersect_fine.orf_start, intersect_fine.orf_end)
    intersect_fine['overlap_end'] = np.where(intersect_fine.strand == '+', intersect_fine.orf_end, intersect_fine.orf_start)
    intersect_fine.drop(columns=['orf_start', 'orf_end'], inplace=True)

    ## total length overlapping is equal to ORF length: overlap len = ORF len
    intersect_fine = pd.merge(intersect_fine, orfs_pred[['orf_id', 'orf_len']], how='inner', on='orf_id')
    intersect_fine = intersect_fine[intersect_fine.overlap_tot == intersect_fine.orf_len]
    ## some ORFs predicted by ribotricer were excluded at this stop because no stop codon (cds_end_nf)

    ## total length overlapping is equal to ORF length: overlap tiv span = ORF len
    intersect_fine['overlap_tstart'] = intersect_fine.apply(lambda r: gtf[r.tx_name].gpos_to_tpos(r.overlap_start)[0], axis=1)
    intersect_fine['overlap_tend'] = intersect_fine.apply(lambda r: gtf[r.tx_name].gpos_to_tpos(r.overlap_end)[0], axis=1)
    intersect_fine['overlap_txlen'] = np.abs(intersect_fine.overlap_tend - intersect_fine.overlap_tstart) + 1
    intersect_fine = intersect_fine[intersect_fine.overlap_txlen == intersect_fine.orf_len]

    # remove ORFs overlapping CDS ends ========================================
    if exclude_cds == True:
        print('...remove ORFs overlapping CDS ends', file=sys.stderr)
        # remove ones that overlap with CDS starts
        intersect_fine = intersect_fine.merge(
            cds_ends.loc[cds_ends.type == 'cds_start', ['chrom', 'gposn']],
            how='outer', indicator=True, left_on=['chrom', 'overlap_start'], right_on=['chrom', 'gposn'])
        intersect_fine = intersect_fine[intersect_fine._merge == 'left_only']
        intersect_fine.drop(columns=['gposn', '_merge'], inplace=True)
        # remove ones that overlap with CDS ends
        intersect_fine = intersect_fine.merge(
            cds_ends.loc[cds_ends.type == 'cds_end', ['chrom', 'gposn']],
            how='outer', indicator=True, left_on=['chrom', 'overlap_end'], right_on=['chrom', 'gposn'])
        intersect_fine = intersect_fine[intersect_fine._merge == 'left_only']
        intersect_fine.drop(columns=['gposn', '_merge'], inplace=True)

    # get complete ncORF candidates ===============================================
    ncorf = intersect_fine[['tx_name', 'orf_id', 'overlap_tstart', 'overlap_tend', 'overlap_start', 'overlap_end', 'orf_len']]
    ncorf.columns = ['tx_name', 'orf_id', 'tstart', 'tend', 'gstart', 'gend', 'orf_len']
    ncorf = ncorf.astype({'tstart': int, 'tend': int, 'gstart': int, 'gend': int, 'orf_len': int})
    # merge with clean transcripts obtained in previous steps
    ncorf = pd.merge(txinfo, ncorf, on='tx_name', how='inner')
    ncorf = ncorf.fillna({'cds_start': 0, 'cds_end': 0}).astype({'cds_start': int, 'cds_end': int})

    # reannotate ORF types ========================================================
    print('...classify ORFs', file=sys.stderr)

    # append simplified transcrpt type
    biotype_table = get_biotype_table()
    ncorf = pd.merge(ncorf, biotype_table, left_on='transcript_biotype', right_on='biotype', how='left')
    ncorf.loc[ncorf.txtype.isna(), 'txtype'] = 'misc'  # in case txtypes are not listed in the biotype table
    ncorf.drop(columns=['biotype'], inplace=True)

    # classification
    ncorf['phase_start'] = (ncorf.tstart - ncorf.cds_start) % 3
    ncorf['phase_end'] = (ncorf.tend - ncorf.cds_end) % 3
    ncorf['orf_type'] = ncorf.apply(orf_typing, axis=1)

    # post-classification filtering and sorting
    if exclude_cds == True:
        misc_cds_variants = ncorf[ncorf.orf_type.isin(
            ['iCDS', 'sCDS', 'N_extension', 'N_truncation', 'C_extension', 'C_truncation', 'wCDS'])]
        ncorf = ncorf[~ncorf.orf_id.isin(misc_cds_variants.orf_id)]

    # sort by orf_id and txtype
    cat_txtype = pd.api.types.CategoricalDtype(
        ['mRNA', 'lncRNA', 'ncRNA', 'pseudogene', 'IGTR', 'misc'], ordered=True)
    ncorf = ncorf.astype({'txtype': cat_txtype})
    ncorf = ncorf.sort_values(by=['orf_id', 'txtype'])

    # output ======================================================================
    print('...write output', file=sys.stderr)

    if pred_method == 'bed12':
        # bed12 has no information about start codon sequences and original ORF types.
        # start codon sequences should be determined manually afterwards
        columns_keep = ['orf_id', 'tx_name', 'txtype', 'orf_type',
                        'tstart', 'tend', 'chrom', 'strand', 'gstart', 'gend', 'orf_len', 'cds_start',
                        'cds_end', 'gene_id', 'gene_name', 'gene_biotype', 'transcript_biotype', 'tx_len',
                        'cds_len', 'utr5_len', 'utr3_len', 'tx_complete', 'phase_start', 'phase_end']
        final_table = ncorf[columns_keep]
    else:
        final_table = pd.merge(ncorf, orfs_pred[['orf_id', 'start_codon']], on='orf_id', how='left')
        final_table = pd.merge(final_table, orfs_pred[['orf_id', 'tx_name', 'orf_type_ori']],
                            on=['orf_id', 'tx_name'], how='left')
        columns_keep = ['orf_id', 'tx_name', 'txtype', 'start_codon', 'orf_type', 'orf_type_ori',
                        'tstart', 'tend', 'chrom', 'strand', 'gstart', 'gend', 'orf_len', 'cds_start',
                        'cds_end', 'gene_id', 'gene_name', 'gene_biotype', 'transcript_biotype', 'tx_len',
                        'cds_len', 'utr5_len', 'utr3_len', 'tx_complete', 'phase_start', 'phase_end']
        final_table = final_table[columns_keep]

    print(final_table.groupby(by='orf_type').size(), file=sys.stderr)
    final_table.to_csv(path_output, sep='\t', index=False)
    return


if __name__ == '__main__':
    main()

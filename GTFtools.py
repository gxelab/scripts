# -*- coding: utf-8 -*-
"""
Operations on genomic intervals stored in GTF file

note:
- all the exons of a gene should be on the same strand. Genes with exons trans-
  spliced from the other strand, like mod(mdg4) in D. melanogaster, should be
  excluded (before or after).
- stop codon is not part of the CDS, at least for Ensembl GTF
-----------------------------------------------------------------
@author: zh (mt1022)
@date: Fri Aug 27 2021
"""
import sys
import re
import gzip
import fileinput
import csv
from dataclasses import dataclass
from typing import List


@dataclass
class Region:
    start: int
    end: int
    
    def __post_init__(self):
        if self.start > self.end:
            raise ValueError('Invalid region boundary!')

    def __len__(self):
        return self.end - self.start + 1


@dataclass
class Gene:
    gene_id: str
    chrom: str = ''
    strand: str = '+'

class Transcript:
    def __init__(self, tx_id: str, gene: Gene):
        self.tx_id: str = tx_id
        self.gene: Gene = gene
        self.exons: List[Region] = []
        self.cdss: List[Region] = []
        self.stop_codon: List[Region] = []
    
    def add_region(self, region, region_type):
        if region_type == 'exon':
            self.exons.append(region)
        elif region_type == 'CDS':
            self.cdss.append(region)
        elif region_type == 'stop_codon':
            self.stop_codon.append(region)
        return
    
    def update(self):
        """
        Update the order of regions so that operations related to intervals can
        work correctly.
        """
        self.exons = sorted(self.exons, key=lambda r: r.start)
        self.cdss = sorted(self.cdss, key=lambda r: r.start)
        self.stop_codon = sorted(self.stop_codon, key=lambda r: r.start)
        return
    
    def __len__(self):
        return sum(len(i) for i in self.exons)
    
    @property
    def n_exons(self):
        return len(self.exons)

    @property
    def cds_len(self):
        return sum(len(i) for i in self.cdss)

    @property
    def tx_start(self):
        if self.gene.strand == '+':
            return self.exons[0].start
        else:
            return self.exons[-1].end

    @property
    def tx_end(self):
        if self.gene.strand == '+':
            return self.exons[-1].end
        else:
            return self.exons[0].start
    
    @property
    def cds_start(self):
        if len(self.cdss) == 0:
            return None
        elif self.gene.strand == '+':
            return self.cdss[0].start
        else:
            return self.cdss[-1].end
    
    @property
    def cds_end(self):
        if len(self.cdss) == 0:
            return None
        elif self.gene.strand == '+':
            return self.cdss[-1].end
        else:
            return self.cdss[0].start
    
    @property
    def stop_codon_start(self):
        if len(self.stop_codon) == 0:
            return None
        elif self.gene.strand == '+':
            return self.stop_codon[-1].end
        else:
            return self.stop_codon[0].start
    
    @property
    def stop_codon_end(self):
        if len(self.stop_codon) == 0:
            return None
        elif self.gene.strand == '+':
            return self.stop_codon[-1].end
        else:
            return self.stop_codon[0].start
    
    @property
    def introns(self):
        if len(self.exons) == 1:
            return []
        else:
            return [Region(self.exons[i].end + 1, self.exons[i+1].start - 1)
                    for i in range(self.n_exons - 1)]

    def tpos_to_gpos(self, pos):
        if pos < 1:
            return 0
        elif pos > len(self):
            return -1
        else:
            if self.gene.strand == '-':
                pos = len(self) - pos + 1
            for i in range(self.n_exons):
                if len(self.exons[i]) < pos:
                    pos -= len(self.exons[i])
                else:
                    return self.exons[i].start + pos - 1

    def gpos_to_tpos(self, pos):
        if pos < self.exons[0].start:
            tpos = self.exons[0].start - pos
            ptype = 'upstream' if self.gene.strand == '+' else 'downstream'
            return tpos, ptype
        elif pos > self.exons[-1].end:
            tpos = pos - self.exons[-1].end
            ptype = 'downstream' if self.gene.strand == '+' else 'upstream'
            return tpos, ptype
        else:
            tpos = 0
            for i in range(self.n_exons):
                if self.exons[i].start <= pos:
                    if self.exons[i].end <= pos:
                        tpos += len(self.exons[i])
                    else:
                        tpos += pos - self.exons[i].start + 1
                else:
                    if self.exons[i-1].end < pos:
                        if self.gene.strand == '+':
                            ptype = 'intron_' + str(i)
                            tpos = pos - self.exons[i - 1].end
                        else:
                            ptype = 'intron_' + str(len(self.exons) - i)
                            tpos = self.exons[i].start - pos
                        return tpos, ptype
                    break
            ptype = 'exon'
            tpos = tpos if self.gene.strand == '+' else len(self) - tpos + 1
            return tpos, ptype
    
    def cpos_to_gpos(self, pos):
        """
        transform CDS coordinate to genomic coordinate
        """
        tpos = self.gpos_to_tpos(self.cds_start)[0] + pos - 1
        gpos = self.tpos_to_gpos(tpos)
        return gpos

    def gpos_to_cpos(self, pos):
        """
        transform genomic coordinate to CDS coordinate
        """
        tpos = self.gpos_to_tpos(pos)[0]
        cpos = tpos - self.gpos_to_tpos(self.cds_start)[0] + 1
        return cpos


    def tiv_to_giv(self, pos1, pos2):
        """
        given transcript region boundary:
        return one or more(for features spanning more than one exon)
        exonic region interval(s) in list of string interval
        """
        cod1 = self.tpos_to_gpos(pos1)
        cod2 = self.tpos_to_gpos(pos2)
        start = min(cod1, cod2)
        end = max(cod1, cod2)
        givs = []

        for i in range(self.n_exons):
            if self.exons[i].end < start:
                continue
            if self.exons[i].start > end:
                break
            if self.exons[i].start <= start:
                if self.exons[i].end <= end:
                    givs.append(Region(start, self.exons[i].end))
                else:
                    givs.append(Region(start, end))
            else:
                if self.exons[i].end <= end:
                    givs.append(Region(self.exons[i].start, self.exons[i].end))
                else:
                    givs.append(Region(self.exons[i].start, end))
        return givs
    
    @property
    def five_prime_utrs(self):
        if len(self.cdss) == 0 or self.cds_start == self.tx_start:
            return []
        else:
            return self.tiv_to_giv(1, self.gpos_to_tpos(self.cds_start)[0] - 1)
    
    @property
    def three_prime_utrs(self):
        if len(self.cdss) == 0 or self.stop_codon_end == self.tx_end or self.cds_end == self.tx_end:
            return []
        else:
            if len(self.stop_codon) > 0:
                return self.tiv_to_giv(self.gpos_to_tpos(self.stop_codon_end)[0] + 1, len(self))
            else:
                return self.tiv_to_giv(self.gpos_to_tpos(self.cds_end)[0] + 1, len(self))
    
    def get_txinfo(self):
        """
        Transcript information
        """
        pass
    
    def format_region_bed(self, r):
        pass
    
    def format_region_bed12(self, rs):
        rs = sorted(rs, key=lambda r: r.start)
        starts = [r.start - 1 for r in rs]
        ends = [r.end for r in rs]

        blockstart = [str(x - starts[0]) for x in starts]
        blocksize = [str(len(r)) for r in rs]

        s = [self.gene.chrom, starts[0], ends[-1], self.tx_id, self.gene.gene_id,
             self.gene.strand, '0', '0', '0', len(starts)]
        s = s + [','.join(blocksize) + ',', ','.join(blockstart) + ',']
        return s

###############################################################################
# functions ###################################################################
###############################################################################
def parse_gtf(gtf_file):
    """
    read GTF file
    """
    gtf = {}
    if gtf_file.endswith('.gz'):
        f = gzip.open(gtf_file, 'rt')
    else:
        f = open(gtf_file)

    for line in f:
        if line[0] == '#':
            continue
        ary = line.strip().split('\t')
        m = re.search(r'gene_id "(.*?)".*?transcript_id "(.*?)"', ary[8])
        if m:
            if m.group(2) in gtf:
                gtf[m.group(2)].add_region(region = Region(int(ary[3]), int(ary[4])), region_type=ary[2])
            else:
                gene = Gene(gene_id=m.group(1), chrom=ary[0], strand=ary[6])
                tx = Transcript(tx_id=m.group(2), gene=gene)
                tx.add_region(region = Region(int(ary[3]), int(ary[4])), region_type=ary[2])
                gtf[m.group(2)] = tx
    f.close()

    for tx in gtf:
        gtf[tx].update()
    return gtf



def exon_to_bed(gtf_file):
    gtf = parse_gtf(gtf_file)
    for tx in gtf:
        items = tx.format_region_bed12(self.exons)
        print('\t'.join(str(i) for i in items))
    return


def cds_to_bed(gtf_file):
    gtf = parse_gtf(gtf_file)
    for tx in gtf:
        if len(tx.cdss) > 0:
            items = tx.format_region_bed12(self.cdss)
            print('\t'.join(str(i) for i in items))
    return


def utr5_to_bed(gtf_file):
    gtf = parse_gtf(gtf_file)
    for tx in gtf:
        tx_utr5 = tx.five_prime_utrs
        if len(tx_utr5) > 0:
            items = tx.format_region_bed12(tx_utr5)
            print('\t'.join(str(i) for i in items))
    return


def utr3_to_bed(gtf_file):
    gtf = parse_gtf(gtf_file)
    for tx in gtf:
        tx_utr3 = tx.three_prime_utrs
        if len(tx_utr3) > 0:
            items = tx.format_region_bed12(tx_utr3)
            print('\t'.join(str(i) for i in items))
    return


def t2g(gtf_file, tfile):
    """
    param tfile: tab-delimited file, 1st column=tx, 2nd column = tpos
    """
    gtf = parse_gtf(gtf_file)
    with open(tfile) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            try:
                tx = gtf[row[0]]
                gpos = tx.tpos_to_gpos(int(row[1]))
                row += [tx.gene.chrom, tx.gene.strand, str(gpos)]
            except KeyError:
                print('Tx isoform {} was not found in GTF file!'.format(row[0]), file=sys.stderr)
                row += ['NA'] * 3
            print('\t'.join(row))
    return


def g2t(gtf_file, gfile):
    """
    param gfile: tab-delimited file, 1st column=tx, 2nd column = gpos
    """
    gtf = parse_gtf(gtf_file)
    with open(gfile) as fh:
        for row in csv.reader(fh, delimiter='\t'):
            try:
                tx = gtf[row[0]]
                tpos, ptype = tx.gpos_to_tpos(int(row[1]))
                row += [str(tpos), ptype]
            except KeyError:
                print('Tx isoform {} was not found in GTF file!'.format(row[0]), file=sys.stderr)
                row += ['NA'] * 2
            print('\t'.join(row))
    return


def tiv2giv(gtf_file, tivfile):
    """
    param tivfile: tab-delimited, first three columns are tx_id, start, and end, 1-based
    """
    gtf = parse_gtf(gtf_file)
    with open(tivfile) as fh:
        for row in csv.reader(fh, delimiter='\t'):
            try:
                tx = gtf[row[0]]
                givs = tx.tiv_to_giv(int(row[1]), int(row[2]))
                print('\t'.join(str(i) for i in tx.format_region_bed12(givs)))
            except KeyError:
                print('Tx isoform {} was not found in GTF file!'.format(row[0]), file=sys.stderr)
    return


def giv2tiv(gtf_file, givfile):
    """
    param givfile: tab-delimited, first three columns are tx_id, start, and end, 1-based
    """
    gtf = parse_gtf(gtf_file)
    with open(givfile) as fh:
        for row in csv.reader(fh, delimiter='\t'):
            try:
                tx = gtf[row[0]]
                if tx.gene.strand == '+':
                    tiv_l = list(tx.gpos_to_tpos(int(row[1])))
                    tiv_r = list(tx.gpos_to_tpos(int(row[2])))
                else:
                    tiv_l = list(tx.gpos_to_tpos(int(row[2])))
                    tiv_r = list(tx.gpos_to_tpos(int(row[1])))
                tiv = [str(tiv_l[0]), str(tiv_r[0]), tiv_l[1], tiv_r[1]]
                row += tiv
            except KeyError:
                row += ['NA'] * 4
                print('Tx isoform {} was not found in GTF file!'.format(row[0]), file=sys.stderr)
            print('\t'.join(row))
    return


def tx_info(gtf_file):
    """
    note: stop codon is counted for CDS length, so that cds + utr5 + utr3 = transcript length
    """
    gtf = parse_gtf(gtf_file)
    header = ['tx_id', 'gene_id', 'chrom', 'strand', 'len', 'len_cds', 'len_utr5', 'len_utr3']
    print('\t'.join(header))
    for tx in gtf:
        out = [tx.tx_id, tx.gene.gene_id, tx.gene.chrom, tx.gene.strand]
        len_tx = len(tx)
        len_utr5 = sum(len(i) for i in tx.five_prime_utrs)
        len_cds = sum(len(i) for i in tx.cdss) + sum(len(i) for i in tx.stop_codon)
        len_utr3 = len_tx - len_cds - len_utr5
        out += [str(i) for i in [len_tx, len_cds, len_utr5, len_utr3]]
        print('\t'.join(out))
    return


#if __name__ == "__main__":
#    parser = argparse.ArgumentParser(
#        prog='GTFtools.py',
#        description='GTF file manipulation')
#    # run mode
#    parser.add_argument(
#        '-M', '--runMode',
#        help='select a subfunction to run',
#        choices=['longest_transcript', 'tssFlank'])
#    # gtf input
#    parser.add_argument(
#        '-g', '--gtf',
#        type=argparse.FileType('r'),
#        help='input gtf file',
#        default=sys.stdin)
#    parser.add_argument(
#        '-l', '--length',
#        type=int,
#        help='length')
#    parser.add_argument(
#        '-d', '--downstreamLen',
#        type=int,
#        help='length of downstream to extract')
#    args = parser.parse_args()

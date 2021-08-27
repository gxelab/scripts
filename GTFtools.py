# -*- coding: utf-8 -*-
"""
Operations on genomic intervals stored in GTF file
--------------------------------------------------------------
@author: zh (mt1022)
@date: Fri Aug 27 2021
"""

import sys
import re
import gzip
import fileinput


class Transcript:
    def __init__(self, mrna_id, gene_id, mrna_chr, strand, desp):
        self.id = mrna_id
        self.gene = gene_id
        self.chr = mrna_chr
        self.strand = strand
        self.exon_starts = []
        self.exon_ends = []
        self.cds_starts = []
        self.cds_ends = []
        self.utr5_starts = []
        self.utr5_ends = []
        self.utr3_starts = []
        self.utr3_ends = []
        self.start_codon = []
        self.stop_codon = []
        self.desp = desp

    def __len__(self):
        length = 0
        nregion = len(self.exon_starts)
        for i in range(nregion):
            length += self.exon_ends[i] - self.exon_starts[i] + 1
        return(length)

    def cds_len(self):
        length = 0
        nregion = len(self.cds_starts)
        for i in range(nregion):
            length += self.cds_ends[i] - self.cds_starts[i] + 1
        return(length)

    def add_feature(self, feature):
        if feature[0] == 'exon':
            self.exon_starts.append(feature[1])
            self.exon_ends.append(feature[2])
        elif feature[0] == 'CDS':
            self.cds_starts.append(feature[1])
            self.cds_ends.append(feature[2])
        elif feature[0] == '5UTR':
            self.utr5_starts.append(feature[1])
            self.utr5_ends.append(feature[2])
        elif feature[0] == '3UTR':
            self.utr3_starts.append(feature[1])
            self.utr3_ends.append(feature[2])
        elif feature[0] == 'start_codon':
            self.start_codon.append([x for x in feature[1:]])
        elif feature[0] == 'stop_codon':
            self.stop_codon.append([x for x in feature[1:]])
        else:
            pass
            # ensembl GTF have 'transcript' and 'gene' feature
            # print("error")
        return

    def get_transcript_start(self):
        """
        get genomic coordinate of transcript start
        """
        pos = -1
        if self.strand == '+':
            pos = sorted(self.exon_starts)[0]
        else:
            pos = sorted(self.exon_ends)[-1]
        return(pos)

    def get_transcript_end(self):
        """
        get genomic coordinate of transcript end
        """
        pos = -1
        if self.strand == '+':
            pos = sorted(self.exon_ends)[-1]
        else:
            pos = sorted(self.exon_starts)[0]
        return(pos)

    def get_cds_start(self):
        """
        get genomic coordinate of CDS start
        """
        cds_starts = self.cds_starts
        cds_ends = self.cds_ends
        cds_starts.sort()
        cds_ends.sort()

        if self.strand == '+':
            return(cds_starts[0])
        else:
            return(cds_ends[-1])

    def get_cds_end(self):
        """
        get genomic coordinate of CDS end
        """
        cds_starts = self.cds_starts
        cds_ends = self.cds_ends
        cds_starts.sort()
        cds_ends.sort()

        if self.strand == '+':
            return(cds_ends[-1])
        else:
            return(cds_starts[0])

    def get_intron_starts(self):
        """
        get start position of all introns
        """
        if self.strand == '+':
            return([int(x) + 1 for x in self.exon_ends[:-1]])
        else:
            return([int(x) - 1 for x in self.exon_starts[1:]])

    def get_intron_ends(self):
        """
        get end position of all introns
        """
        if self.strand == '+':
            return([int(x) - 1 for x in self.exon_starts[1:]])
        else:
            return([int(x) + 1 for x in self.exon_ends[:-1]])

    def tpos_to_gpos(self, pos):
        """
        transform transcript coordinate to genomic coordinate
        """
        index = pos
        nregion = len(self.exon_starts)
        exon_starts = sorted(self.exon_starts)
        exon_ends = sorted(self.exon_ends)
        if self.strand == '+':
            for i in range(nregion):
                if exon_ends[i] - exon_starts[i] + 1 < index:
                    index -= exon_ends[i] - exon_starts[i] + 1
                elif exon_ends[i] - exon_starts[i] + 1 >= index:
                    return(exon_starts[i] + index - 1)
        else:
            for i in reversed(range(nregion)):
                if exon_ends[i] - exon_starts[i] + 1 < index:
                    index -= exon_ends[i] - exon_starts[i] + 1
                elif exon_ends[i] - exon_starts[i] + 1 >= index:
                    return(exon_ends[i] - index + 1)

    def gpos_to_tpos(self, pos):
        """
        transform genomic coordinate to transcript coordinate
        """
        dist = 0
        nregion = len(self.exon_starts)
        if self.strand == '+':
            for i in range(nregion):
                if self.exon_starts[i] <= pos <= self.exon_ends[i]:
                    dist += pos - self.exon_starts[i] + 1
                elif self.exon_ends[i] < pos:
                    dist += self.exon_ends[i] - self.exon_starts[i] + 1
            return(dist)
        else:
            for i in range(nregion):
                if self.exon_starts[i] <= pos <= self.exon_ends[i]:
                    dist += self.exon_ends[i] - pos + 1
                elif self.exon_starts[i] > pos:
                    dist += self.exon_ends[i] - self.exon_starts[i] + 1
            return(dist)

    def cpos_to_gpos(self, pos):
        """
        transform CDS coordinate to genomic coordinate
        """
        tpos = self.gpos_to_tpos(self.get_cds_start()) + pos - 1
        gpos = self.tpos_to_gpos(tpos)
        return(gpos)

    def gpos_to_cpos(self, pos):
        """
        transform genomic coordinate to CDS coordinate
        """
        tpos = self.gpos_to_tpos(pos)
        cpos = tpos - self.gpos_to_tpos(self.get_cds_start()) + 1
        return(cpos)

    def print_sregion(self, pos1, pos2):
        """
        return simple region in string interval
        """
        s = [str(pos1), str(pos2)]
        return(s)

    def print_cregion(self, pos1, pos2):
        """
        given transcript region boundary:
        return one or more(for features spanning more than one exon)
        exonic region interval(s) in list of string interval
        """
        cod1 = self.tpos_to_gpos(pos1)
        cod2 = self.tpos_to_gpos(pos2)
        start = min(cod1, cod2)
        end = max(cod1, cod2)
        nregion = len(self.exon_starts)
        res = []

        # print(pos1, pos2, 'start: ', start, 'end', end)
        for i in range(nregion):
            if start >= self.exon_starts[i] and end <= self.exon_ends[i]:
                left = start
                right = end
                res.append(self.print_sregion(left, right))
            elif start < self.exon_starts[i] and self.exon_starts[i] <= end <= self.exon_ends[i]:
                left = self.exon_starts[i]
                right = end
                res.append(self.print_sregion(left, right))
            elif end > self.exon_ends[i] and self.exon_starts[i] <= start <= self.exon_ends[i]:
                left = start
                right = self.exon_ends[i]
                res.append(self.print_sregion(left, right))
            elif start <= self.exon_starts[i] and end >= self.exon_ends[i]:
                left = self.exon_starts[i]
                right = self.exon_ends[i]
                res.append(self.print_sregion(left, right))
        return(res)

    def is_valid_tpos(self, pos):
        """
        is input a valid transcript coordinate for this mRNA
        """
        if pos < 1 or pos > len(self):
            return(False)
        else:
            return(True)

    def is_valid_gpos(self, pos):
        """
        is input a valid genomic coordinate for this mRNA
        """
        exon_starts = self.exon_starts
        exon_ends = self.exon_ends
        exon_starts.sort()
        exon_ends.sort()
        if pos < exon_starts[0] or pos > exon_ends[-1]:
            return(False)
        else:
            return(True)

    def is_coding(self):
        """
        is this mRNA protein-coding?
        """
        if len(self.cds_starts) == 0:
            return(False)
        else:
            return(True)


#==============================================================================
# functions
#==============================================================================


def GTFparse(gtf_file):
    """
    read GTF file
    """
    all_transcript = {}
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
            if m.group(2) in all_transcript:
                all_transcript[m.group(2)].add_feature([ary[2], int(ary[3]), int(ary[4])])
            else:
                mrna = Transcript(m.group(2), m.group(1), ary[0], ary[6], ary[8])
                mrna.add_feature([ary[2], int(ary[3]), int(ary[4])])
                all_transcript[m.group(2)] = mrna
    f.close()
    return(all_transcript)


def longest_transcript(gtf_file, p=True):
    """
    print longest transcript info for each gene

    longest_transcript(str gtf_file, logic p) -> dict {gene:longest_mRNA}
    """
    gtf = GTFparse(gtf_file)
    all_gene = {}
    # read transcripts length for each gene
    for mrna in gtf:
        if gtf[mrna].gene in all_gene:
            all_gene[gtf[mrna].gene][mrna] = len(gtf[mrna])
        else:
            all_gene[gtf[mrna].gene] = {mrna: len(gtf[mrna])}
    # get longest transcript for each gene
    longest_mrna = {}
    all_mrna = {}  # hash for fast search in printing step
    for gene in all_gene:
        longest = ''
        longest_len = 0
        for mrna in all_gene[gene]:
            all_mrna[mrna] = 0
            if all_gene[gene][mrna] > longest_len:
                if longest != '':
                    all_mrna[longest] = 0
                longest = mrna
                longest_len = all_gene[gene][mrna]
                all_mrna[longest] = 1
        #print(gene + '\t' + longest)
        longest_mrna[gene] = longest

    # print annotation of longest transcript
    set_mrna = set(longest_mrna.values())
    set_gene = set(longest_mrna.keys())
    if p:
        for line in open(gtf_file):
            if line[0] == '#':
                continue
            ary = line.strip().split('\t')
            if ary[2] == 'gene':  # if ary[2] == 'gene'
                m = re.search(r'gene_id "(.*?)"', ary[8])
                if m and m.group(1) in set_gene:
                    print(line, end='')
                continue
            m = re.search(r'gene_id "(.*?)".*?transcript_id "(.*?)"', ary[8])
            if m and m.group(2) in set_mrna:  # if all_mrna[m.group(2)]:
                print(line, end='')
    return(longest_mrna)


def tss_downstream(gtf_file, nt):
    """
    get annotation of tss downstream nt region for longest mRNA of each gene

    tss_downstream(str, list, int) -> list

    status: manually confirmed on 2015-09-18
    """
    longest_mrna = longest_transcript(gtf_file, p=False)
    gtf = GTFparse(gtf_file)
    for geneID in longest_mrna:
        mrnaID = longest_mrna[geneID]
        mrna = gtf[mrnaID]
        if mrna.strand == '+':
            rstart = mrna.get_transcript_start()
            rend = rstart + nt - 1
        else:
            rend = mrna.get_transcript_start()
            rstart = rend - (nt - 1)
        s = [mrna.chr, 'GTFtools', 'tssFlank', str(rstart), str(rend), '.',
             mrna.strand, '.', mrna.desp]
        print("\t".join(s))
    return(0)


def utr5_bed12(gtf_file):
    """
    print UTR5 annotation in bedtools.bed12 format
    """
    gtf = GTFparse(gtf_file)
    for mrnaID in gtf:
        mrna = gtf[mrnaID]
        #print(mrnaID, mrna.gene)
        if len(mrna.cds_starts) == 0:
            continue
        cds_start = mrna.get_cds_start()
        utr5_len = mrna.gpos_to_tpos(cds_start) - 1
        if utr5_len < 1:
            continue
        utr5_interval = mrna.print_cregion(1, utr5_len)
        starts = sorted([int(x[0]) - 1 for x in utr5_interval])
        ends = sorted([int(x[1]) for x in utr5_interval])

        blockstart = [str(x - starts[0]) for x in starts]
        blocksize = [str(ends[x] - starts[x]) for x in range(len(starts))]

        s = [mrna.chr, starts[0], ends[-1], mrna.id, mrna.gene, mrna.strand,
             '0', '0', '0', len(starts)]
        s = s + [','.join(blocksize) + ',', ','.join(blockstart) + ',']
        s = [str(x) for x in s]
        print('\t'.join(s))
    return


def utr5_bed12_ex(gtf_file, up=15, down=15):
    """
    print UTR5 annotation in bedtools.bed12 format with flanking region
    up: bases upstream of TSS in genomic sequence
    down: bases downstream of 5' UTR in CDS sequence
    """
    gtf = GTFparse(gtf_file)
    up = int(up)
    down = int(down)
    for mrnaID in gtf:
        mrna = gtf[mrnaID]
        #print(mrnaID, mrna.gene)
        if len(mrna.cds_starts) == 0:
            continue
        cds_start = mrna.get_cds_start()
        utr5_len = mrna.gpos_to_tpos(cds_start) - 1
        if utr5_len < 1 or mrna.cds_len() < down:
            continue
        # extend to downstream
        utr5_interval = mrna.print_cregion(1, utr5_len + down)
        starts = sorted([int(x[0]) - 1 for x in utr5_interval])
        ends = sorted([int(x[1]) for x in utr5_interval])
        # extends to upstream
        if mrna.strand == '+':
            starts[0] = starts[0] - up
        else:
            ends[-1] = ends[-1] + up
        blockstart = [str(x - starts[0]) for x in starts]
        blocksize = [str(ends[x] - starts[x]) for x in range(len(starts))]

        s = [mrna.chr, starts[0], ends[-1], mrna.id, mrna.gene, mrna.strand,
             '0', '0', '0', len(starts)]
        s = s + [','.join(blocksize) + ',', ','.join(blockstart) + ',']
        s = [str(x) for x in s]
        print('\t'.join(s))
    return


def utr3_bed12(gtf_file):
    """
    print UTR3 annotation in bedtools.bed12 format
    """
    gtf = GTFparse(gtf_file)
    for mrnaID in gtf:
        mrna = gtf[mrnaID]
        #print(mrnaID, mrna.gene)
        if len(mrna.cds_starts) == 0:
            continue
        # in ens gtf: stop codon is not part of CDS
        cds_end = mrna.gpos_to_tpos(mrna.get_cds_end()) + 3
        utr3_len = len(mrna) - cds_end
        if utr3_len < 1:
            continue
        utr3_interval = mrna.print_cregion(cds_end + 1, len(mrna))
        starts = sorted([int(x[0]) - 1 for x in utr3_interval])
        ends = sorted([int(x[1]) for x in utr3_interval])

        blockstart = [str(x - starts[0]) for x in starts]
        blocksize = [str(ends[x] - starts[x]) for x in range(len(starts))]

        s = [mrna.chr, starts[0], ends[-1], mrna.id, mrna.gene, mrna.strand,
             '0', '0', '0', len(starts)]
        s = s + [','.join(blocksize) + ',', ','.join(blockstart) + ',']
        s = [str(x) for x in s]
        print('\t'.join(s))
    return


def cds_bed12(gtf_file):
    """
    print CDS annotation in bedtools.bed12 format

    status: manually confirmed on 2016-0215
    """
    gtf = GTFparse(gtf_file)
    for i in gtf:
        record = gtf[i]
        if len(record.cds_starts) == 0:
            continue
        starts = sorted([x - 1 for x in record.cds_starts])
        ends = sorted(record.cds_ends)
        blockstart = [str(x - starts[0]) for x in starts]
        blocksize = [str(ends[x] - starts[x]) for x in range(len(starts))]

        s = [record.chr, starts[0], ends[-1], record.id, record.gene, record.strand,
             '0', '0', '0', len(starts)]
        s = s + [','.join(blocksize) + ',', ','.join(blockstart) + ',']
        s = [str(x) for x in s]
        print('\t'.join(s))
    return


def cds_bed12_crop(gtf_file, l=0, r=0):
    """
    print truncated (arbitrary length from both ends) CDS annotation

    status: manually confirmed on 2016-0215
    """
    l = int(l)
    r = int(r)
    gtf = GTFparse(gtf_file)
    for i in gtf:
        record = gtf[i]
        if len(record.cds_starts) == 0:
            continue
        t_start = record.gpos_to_tpos(record.get_cds_start()) + l
        t_end = record.gpos_to_tpos(record.get_cds_end()) - r
        #pdb.set_trace()
        if t_start < t_end:
            crop_intervals = record.print_cregion(t_start, t_end)
            starts = sorted([int(x[0]) - 1 for x in crop_intervals])
            ends = sorted([int(x[1]) for x in crop_intervals])
            blockstart = [str(x - starts[0]) for x in starts]
            blocksize = [str(ends[x] - starts[x]) for x in range(len(starts))]
            s = [record.chr, starts[0], ends[-1], record.id, record.gene, record.strand,
                 '0', '0', '0', len(starts)]
            s = s + [','.join(blocksize) + ',', ','.join(blockstart) + ',']
            s = [str(x) for x in s]
            print('\t'.join(s))
    return


def tiv_to_giv(gtf_file, tiv_file):
    """
    convert transcript interval to genomic interval
    :param gtf_file: gtf
    :param tiv_file: tab delimited file (isoform, start, end), both start & end are 1-based;
    :return: bed12
    """
    gtf = GTFparse(gtf_file)
    for line in open(tiv_file, 'r'):
        ary = line.rstrip().split('\t')
        isoform = ary[0]
        t_start = int(ary[1])
        t_end = int(ary[2])
        if isoform not in gtf:
            continue
        record = gtf[isoform]
        if t_start < t_end <= len(record):
            crop_intervals = record.print_cregion(t_start, t_end)
            starts = sorted([int(x[0]) - 1 for x in crop_intervals])
            ends = sorted([int(x[1]) for x in crop_intervals])
            blockstart = [str(x - starts[0]) for x in starts]
            blocksize = [str(ends[x] - starts[x]) for x in range(len(starts))]
            s = [record.chr, starts[0], ends[-1], record.id, record.gene, record.strand,
                 '0', '0', '0', len(starts)]
            s = s + [','.join(blocksize) + ',', ','.join(blockstart) + ',']
            s = [str(x) for x in s]
            print('\t'.join(s))
    return


def tiv_to_giv2(gtf_file, tiv_file):
    """
    convert transcript interval to genomic interval
    :param gtf_file: gtf
    :param tiv_file: tab delimited file (isoform, start, end), both start & end are 1-based;
    :return: bed12
    """
    gtf = GTFparse(gtf_file)
    for line in fileinput.input(tiv_file):
        ary = line.rstrip().split('\t')
        isoform = ary[0]
        t_start = int(ary[1])
        t_end = int(ary[2])
        if isoform not in gtf:
            continue
        record = gtf[isoform]
        if t_start < t_end <= len(record):
            crop_intervals = record.print_cregion(t_start, t_end)
            starts = sorted([int(x[0]) - 1 for x in crop_intervals])
            ends = sorted([int(x[1]) for x in crop_intervals])
            blockstart = [str(x - starts[0]) for x in starts]
            blocksize = [str(ends[x] - starts[x]) for x in range(len(starts))]
            s = [record.chr, starts[0], ends[-1], ary[3], record.id, record.strand,
                 '0', '0', '0', len(starts)]
            s = s + [','.join(blocksize) + ',', ','.join(blockstart) + ',']
            s = [str(x) for x in s]
            print('\t'.join(s))
    return


def convert_g2t(gtf_file, gpos_file):
    """
    convert genomic coordinates to transcript coordinates
    :param gtf_file: gtf
    :param tiv_file: tab delimited file (isoform, gpos:1-based);
    :return: tab-delimited file with an additional column of tpos
    """
    gtf = GTFparse(gtf_file)
    for line in fileinput.input(gpos_file):
        ary = line.rstrip().split('\t')
        isoform = ary[0]
        gpos = int(ary[1])
        if isoform not in gtf:
            tpos = -1
        else:
            record = gtf[isoform]
            tpos = record.gpos_to_tpos(gpos)
        print('\t'.join(ary + [str(tpos)]))
    return


def transcript_bed(gtf_file):
    """
    get bed annotation of transcript start and end and print entries sorted
    by coordinate;
    """
    gtf = GTFparse(gtf_file)
    bed = []
    for i in gtf:
        record = gtf[i]
        if record.strand == '+':
            bed_start = record.get_transcript_start() - 1
            bed_end = record.get_transcript_end()
        else:
            bed_start = record.get_transcript_end() - 1
            bed_end = record.get_transcript_start()
        s = [record.chr, str(bed_start), str(bed_end), i, '0', record.strand]
        bed.append(s)
    bed = sorted(bed, key=lambda x: (x[0], x[1]))
    for i in bed:
        print('\t'.join(i))
    return


def transcript_bed12(gtf_file):
    """
    print transcript annotation in bedtools.bed12 format
    """
    gtf = GTFparse(gtf_file)
    for i in gtf:
        record = gtf[i]
        starts = sorted([x - 1 for x in record.exon_starts])
        ends = sorted(record.exon_ends)
        blockstart = [str(x - starts[0]) for x in starts]
        blocksize = [str(ends[x] - starts[x]) for x in range(len(starts))]

        s = [record.chr, starts[0], ends[-1], record.id, record.gene, record.strand,
             '0', '0', '0', len(starts)]
        s = s + [','.join(blocksize) + ',', ','.join(blockstart) + ',']
        s = [str(x) for x in s]
        print('\t'.join(s))
    return


def gene_info(gtf_file):
    """
    print length information of genes;
    """
    gtf = GTFparse(gtf_file)
    for i in gtf:
        out = [i, gtf[i].gene, str(len(gtf[i]))]
        print('\t'.join(out))
    return

def convert_t2g(gtf_file, coord_file):
    """
    convert transcript coordinates to genomic coordinates
    """
    gtf = GTFparse(gtf_file)
    for line in open(coord_file, 'r'):
        ary = line.rstrip().split()
        isoform = ary[0]
        tpos = ary[1]
        try:
            gpos = gtf[isoform].tpos_to_gpos(int(tpos))
            chrm = gtf[isoform].chr
            strand = gtf[isoform].strand
            gpos = str(gpos)
        except KeyError:
            gpos = 'isoform_missing'
            chrm = 'NA'
            strand = 'NA'
        out = [isoform, tpos, chrm ,strand, gpos]
        print('\t'.join(out))
    return


def convert_t2g_v2(gtf_file, coord_file, col_trans = 1, col_tpos = 2):
    """
    convert transcript coordinates to genomic coordinates
    """
    gtf = GTFparse(gtf_file)
    col_trans = int(col_trans)
    col_tpos = int(col_tpos)
    for line in open(coord_file, 'r'):
        ary = line.rstrip().split()
        isoform = ary[col_trans - 1]
        tpos = ary[col_tpos - 1]
        try:
            chrm = gtf[isoform].chr
            gpos = gtf[isoform].tpos_to_gpos(int(tpos))
            gpos = str(gpos)
            strand = gtf[isoform].strand
        except KeyError:
            chrm = 'chrm'
            gpos = 'gpos'
            strand = 'strand'
        ary += [chrm, gpos, strand]
        print('\t'.join(ary))
    return


def transcript_len_stat(gtf_file):
    """
    summary transcript length, cds length, utr length of each gene in gtf;
    output isoform length, cds length, cds t_start, cds t_end;
    """
    gtf = GTFparse(gtf_file)
    for i in gtf:
        isoform = gtf[i]
        len_trans = len(isoform)
        if isoform.cds_starts != []:
            cds_start = isoform.gpos_to_tpos(isoform.get_cds_start())
            cds_end = isoform.gpos_to_tpos(isoform.get_cds_end())
            len_cds = isoform.cds_len()
            out = [i, isoform.gene, isoform.chr, isoform.strand, len_trans, len_cds, cds_start, cds_end]
        else:
            out = [i, isoform.gene, isoform.chr, isoform.strand, len_trans, 'NA', 'NA', 'NA']
        print('\t'.join(str(j) for j in out))
    return


def convert_t2g_all(gtf_file):
    """
    convert coordinates of all transcript to genomic coordinates
    """
    gtf = GTFparse(gtf_file)
    for isof_id in gtf:
        isoform = gtf[isof_id]
        isof_len = len(isoform)
        if isoform.cds_starts != []:
            cds_start = isoform.gpos_to_tpos(isoform.get_cds_start())
        else:
            cds_start = 'NA'
        for tpos in range(1, isof_len + 1):
            if cds_start != 'NA':
                cpos = tpos - cds_start + 1
            else:
                cpos = 'NA'
            gpos = isoform.tpos_to_gpos(tpos)
            out = [isof_id, str(tpos), str(cpos), isoform.chr, isoform.strand, str(gpos)]
            print('\t'.join(out))
    return


#==============================================================================
# main function
#==============================================================================
if __name__ == '__main__':
    if sys.argv[1] == 'longest_transcript':
        longest_transcript(*sys.argv[2:])
    elif sys.argv[1] == 'utr5_bed12':
        utr5_bed12(*sys.argv[2:])
    elif sys.argv[1] == 'cds_bed12':
        cds_bed12(*sys.argv[2:])
    elif sys.argv[1] == 'utr3_bed12':
        utr3_bed12(*sys.argv[2:])
    elif sys.argv[1] == 'cds_bed12_crop':
        cds_bed12_crop(*sys.argv[2:])
    elif sys.argv[1] == 'transcript_bed':
        transcript_bed(*sys.argv[2:])
    elif sys.argv[1] == 'transcript_bed12':
        transcript_bed12(*sys.argv[2:])
    elif sys.argv[1] == 'utr5_bed12_ex':
        utr5_bed12_ex(*sys.argv[2:])
    elif sys.argv[1] == 'gene_info':
        gene_info(*sys.argv[2:])
    elif sys.argv[1] == 'convert_t2g':
        convert_t2g(*sys.argv[2:])
    elif sys.argv[1] == 'convert_t2g_v2':
        convert_t2g_v2(*sys.argv[2:])
    elif sys.argv[1] == 'convert_g2t':
        convert_g2t(*sys.argv[2:])
    elif sys.argv[1] == 'convert_t2g_all':
        convert_t2g_all(*sys.argv[2:])
    elif sys.argv[1] == 'len_stat':
        transcript_len_stat(*sys.argv[2:])
    elif sys.argv[1] == 'tiv_to_giv':
        tiv_to_giv(*sys.argv[2:])
    elif sys.argv[1] == 'tiv_to_giv2':
        tiv_to_giv2(*sys.argv[2:])
    else:
        print('wrong parameters!')



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

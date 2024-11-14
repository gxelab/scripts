# -*- coding: utf-8 -*-
"""
A python script to infer ORF origination node and mechanism and calculate BLS
score given ancestral sequence reconstruction.
-----------------------------------------------------------------
@author: Lei TY & mt1022
@date: 2024-01-18
"""

from math import floor
from collections import namedtuple
from ete4 import Tree, PhyloTree
from Bio import AlignIO, Seq
import click

def to_dict(aln):
    a2d = dict()
    for record in aln:
        a2d[record.id] = record.seq
    return a2d


def remove_gaps(seq):
    """Remove gaps in a sequence"""
    return Seq.Seq(str(seq).replace('-', ''))


def to_codon(seq):
    codons = []
    for i in range(0, len(seq), 3):
        codons.append(seq[i:i + 3])
    return codons


def find_nth_nonx(string, n, x='-'):
    """
    Find the nth occurrence of a non-x character
    
    Retruns:
      1-based index of nth occurrence if found else -1
    """
    index = 0
    count = 0
    while count < n:
        while string[index] == x:
            index += 1
            if index == len(string):
                return -1
        count += 1
        index += 1
        if index == len(string):
            return -1
    return index


def is_conserved(s1, s2, cutoff=0.7):
    """
    Check whether s2 has conserved ORF structure as s1

    criteria (similar to that in Sandman et al., 2022, MC):
    - start codon at the same position or the next two in-frame codons
    - 70% sequence of the sequence not interrupted by stop codons truncating the ORF
    """
    # whether one of the first codons is a start codon
    tis = ['ATG', 'CTG', 'GTG', 'TTG']
    s2_codons = to_codon(remove_gaps(s2).upper())
    has_start = [codon in tis for codon in s2_codons[:3]]
    if sum(has_start) > 0:
        which_start = has_start.index(True)
    else:
        which_start = 0

    # detemine minimum floor(70% of amino acid length from start of s2)
    aa70 = floor((len(remove_gaps(s1))/3 - 1) * cutoff)
    min_len = (which_start + aa70) * 3
    # find premature stop codons minimum conserved regions
    min_len_start = find_nth_nonx(s2, which_start * 3)
    min_len_stop = find_nth_nonx(s1, min_len)
    s2_orf = remove_gaps(s2[min_len_start:min_len_stop])
    pmsc = s2_orf.translate().find('*') >= 0

    # summarize results
    Result = namedtuple('Result', ['has_start', 'no_pmsc', 'conserved'])
    conserved = (sum(has_start) > 0 and pmsc == 0) + 0
    return Result(sum(has_start) > 0, pmsc == 0, conserved)


def blsum(tree):
    """
    get total branch length of a tree (the root node is ignored)
    """
    return sum([n.dist for n in tree.traverse()][1:])


def lastn(aln, n):
    """index the start (in reverse order) of last n bases"""
    i = 0
    while n > 0 and i < len(aln):
        i += 1
        if aln[-i] != '-':
            n -= 1
    return i


@click.command(context_settings=dict(help_option_names=['-h', '--help'], show_default=True))
@click.argument('path_asr', type=click.STRING)
@click.argument('tree_asr', type=click.STRING)
@click.argument('tree_full', type=click.STRING)
@click.argument('ref_sp', type=click.STRING)
@click.option('-o', '--output_prefix', type=click.STRING, default='_tmp',
              help='output prefix')
@click.option('-c', '--clean', is_flag=True, default=False, help="only write the results file")
def orf_bls(path_asr, tree_asr, tree_full, ref_sp, output_prefix, clean=False):
    """
    Infer age, origination mechanism and Branch Length Score (BLS) of an ORF

    Arguments\n
    path_asr: path to multiple sequence aligment containing ancestral sequence reconstruction (FastML results).\n
    tree_asr: path to tree whose nodes and leaves have associated sequence in alignment (asr).\n
    tree_full: path to global tree with the local_tree as a subtree.\n
    ref_sp: name of the reference species in tree/fasta.\n
    """
    # read ancestral sequence reconstruction results: tree + alignment (including ancestor nodes)
    tasr = PhyloTree(open(tree_asr), parser=1, alignment=path_asr, alg_format='fasta')

    # check whether the ORF of reference species is conserved at each node in the asr tree
    ref_node = tasr[ref_sp]
    for n in tasr.traverse():
        n_orf = is_conserved(ref_node.props['sequence'], n.props['sequence']).conserved
        tasr[n.name].add_prop('orf', n_orf)

    # determine origination node and manner -------------------------------------------------------
    anc_name = [n.name for n in tasr[ref_sp].ancestors()]  # bottom-up list of ancestors of the reference species
    anc_worf = [n.get_prop('orf') for n in tasr[ref_sp].ancestors()]  # whether each ancestor has conserved ORF struct

    # origin_age: the most ancient node among ancestors of the reference species that has conserved ORF structure
    # origin_manner: a node with reconstructed sequence predates the node of origin
    # ref:
    #   Sandmann et al., 2023 MC:
    #  "the most ancient ancestral primatomorpha sequence predicted to contain the full ORF structure"
    #   Vakirlis et al., 2022 CR: 
    #  "the most ancient human ancestor in which at least 70% of the reconstructed ancestral sequence was an intact ORF"
    if anc_worf[-1] == 0:
        origin_manner = 'denovo'
        try:
            origin_age = anc_name[-1 - anc_worf[::-1].index(1)]
        except:
            origin_age = ref_sp
    else:
        origin_manner = 'nondenovo'
        origin_age  = anc_name[-1]

    # get lists of leaves 
    # all leaves of the origin node
    origin_leaves = [n.name for n in tasr[origin_age].leaves()]
    # leaves of the origin node that have conserved ORF structure
    origin_leaves_orf = [n.name for n in tasr[origin_age]  if n.get_prop('orf') == 1]
    # all the leaves that have conserved ORF among species used for ancestral sequence reconstruction (ASR)
    local_leaves_orf = [n.name for n in tasr if n.get_prop('orf') == 1]

    # read the full tree and assign names to the internal node ------------------------------------
    tfull = Tree(open(tree_full), parser=1)
    n = 1
    for node in tfull.traverse('preorder'):
        if not node.is_leaf:
            node.name = f'N{n}'
            n += 1

    # BLS calculation method 1 --------------------------------------------------------------------
    # determine the subtree with the origination node inferred from ASR
    # locate the origination node in the full tree
    origin_age_global = tfull.common_ancestor(origin_leaves).name
    # get the total branch length of the full tree
    bl_all = blsum(tfull)
    # get the total branch length of the subtree tree defined by the origination node
    bl_sub_origin = blsum(tfull[origin_age_global])

    # get the total branch length of connections between nodes with ORFs in the subtree
    # if the is only one species with ORF (i.e., the reference node), bl_orf is the distance
    #   between the origination node and the ref_sp
    # else bl_orf is total branch lengths of the minimal tree connnecting leaves that are descendants
    #   of the origination node, (represented by the MRCA of these leaves), plus the distance between 
    #   the origination node and their MRCA
    if origin_leaves_orf == [ref_sp]:
        bl_orf_origin = tfull.get_distance(origin_age_global, ref_sp)
    else:
        origin_leaves_orf_mrca = tfull.common_ancestor(origin_leaves_orf).name
        tmrca = tfull.copy()
        tmrca.prune(origin_leaves_orf, preserve_branch_length=True)
        bl_orf_origin = (tfull.get_distance(origin_age_global, origin_leaves_orf_mrca) +
                         sum([n.dist for n in tmrca.traverse()][1:]))
    
    # global bls: bl_orf / bl_fulltree
    # local bls: bl_orf / bl_subtree
    bls_global_origin = bl_orf_origin / bl_all
    bls_local_origin = bl_orf_origin / bl_sub_origin if bl_sub_origin != 0 else 0

    # BLS calculation method 2 --------------------------------------------------------------------
    # a naive but intuitive way: determine the origination node as MRCA of all leaves with conserved
    #   ORF structure
    # bl_sub: the total branch length of the subtree defined by naive_age_global
    # bl_orf: total branch length of the subtree defined by leaves with conserved ORF structure， whose
    # root node is the same as naive_age_global
    naive_age_global = tfull.common_ancestor(local_leaves_orf).name
    if local_leaves_orf == [ref_sp]:
        bl_sub_naive = 0
        bl_orf_naive = 0
    else:
        bl_sub_naive = sum([n.dist for n in tfull[naive_age_global].traverse()][1:])
        tnaive = tfull[naive_age_global].copy()
        tnaive.prune(local_leaves_orf, preserve_branch_length=True)
        bl_orf_naive = sum([n.dist for n in tnaive.traverse()][1:])
    
    bls_global_naive = bl_orf_naive / bl_all
    bls_local_naive = bl_orf_naive / bl_sub_naive if bl_sub_naive != 0 else 0

    # save results
    out = [origin_manner, origin_age, origin_age_global, bl_all,
           bl_sub_origin, bl_orf_origin, bls_global_origin, bls_local_origin,
           naive_age_global, bl_sub_naive, bl_orf_naive, bls_global_naive, bls_local_naive]
    field = ['origin_manner', 'origin_age', 'origin_age_global', 'bl_all',
             'bl_sub_origin', 'bl_orf_origin', 'bls_global_origin', 'bls_local_origin',
             'naive_age_global', 'bl_sub_naive', 'bl_orf_naive', 'bls_global_naive', 'bls_local_naive']
    with open(output_prefix + '.result.txt', 'wt') as fh:
        print('\t'.join(f'{i}\t{j}' for i, j in zip(field, out)), file=fh)

    # other information for debugging
    if clean == False:
        with open(output_prefix + '.node_status_plot.txt', 'wt') as fh:
            print('Tree ASR:', end='\n', file=fh) 
            print(tasr.to_str(compact=True, props=['name', 'orf']), file=fh)
            print('Tree Global:', end='\n', file=fh) 
            print(tfull.to_str(compact=True, props=['name']), file=fh)
        # newick for Newick format; nhx for New Hampshire eXtended format
        with open(output_prefix + '.node_status.nhx', 'wt') as fh:
            print(tasr.write(props=['orf']), file=fh)
        # alignment of complete ORFs under origination node
        sv_lastn = lastn(tasr[ref_sp].props['sequence'], 3)
        if origin_leaves_orf != [ref_sp]:
            with open(output_prefix + '.orfs_origin.fa', 'wt') as fh:
                for n in origin_leaves_orf:
                    seq = tasr[n].props['sequence'][:-sv_lastn]
                    print(f'>{n}\n{seq}', end = '\n', file=fh)
            with open(output_prefix + '.orfs_origin.nwk', 'wt') as fh:
                print(tmrca.write(), file=fh)
        # alignment of all complete ORFs
        if local_leaves_orf != [ref_sp]:
            with open(output_prefix + '.orfs_naive.fa', 'wt') as fh:
                for n in local_leaves_orf:
                    seq = tasr[n].props['sequence'][:-sv_lastn]
                    print(f'>{n}\n{seq}', end = '\n', file=fh)
            with open(output_prefix + '.orfs_naive.nwk', 'wt') as fh:
                print(tnaive.write(), file=fh)
    return out


if __name__ == '__main__':
    orf_bls()

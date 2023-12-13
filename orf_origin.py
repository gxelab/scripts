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


@click.command(context_settings=dict(help_option_names=['-h', '--help'], show_default=True))
@click.argument('asr', type=click.STRING)
@click.argument('local_tree', type=click.STRING)
@click.argument('global_tree', type=click.STRING)
@click.option('-o', '--output_prefix', type=click.STRING, default='_tmp')
def infer_orf_age(asr, local_tree, global_tree, output_prefix):
    """
    Infer age and origin mechanism of an ORF

    Arguments
    asr: path to multiple sequence aligment containing ancestral sequence reconstruction (FastML results).
    local_tree: path to tree whose nodes and leaves have associated sequence in alignment (asr).
    global_tree: path to global tree with the local_tree as a subtree.
    output_prefix: output_prefix.
    """
    aln = AlignIO.read(asr, 'fasta')
    tree = Tree(open(local_tree))

    # evaluate whether ORF structure is conserved at each node/leaf
    aln_dict = to_dict(aln)
    for node in tree.traverse():
        node_status = is_conserved(aln_dict['hg38'], aln_dict[node.name])
        tree[node.name].add_prop('orf', node_status.conserved)
    
    # determine the local age and origination mechanism of the human ORF
    ancestors = [node.name for node in tree['hg38'].ancestors()]
    orf_status = [node.props['orf'] for node in tree['hg38'].ancestors()]
    
    try:
        first_non = orf_status.index(0)
        if first_non == 0:
            local_age = 'hg38'
        else:
            local_age = ancestors[first_non - 1]
        orf_origin = 'denovo'
    except:
        local_age = ancestors[-1]
        orf_origin = 'non_denovo'

    # determine the age of ORF in the global tree
    global_tree = Tree(open(global_tree))
    n = 1
    for node in global_tree.traverse(strategy='preorder'):
        if node.name == None:
            node.name = f'N{n}'
            n += 1
    orf_leaves = [node.name for node in tree[local_age].leaves()]
    global_age = global_tree.common_ancestor(orf_leaves).name

    # save results
    with open(output_prefix + '.node_status_plot.txt', 'wt') as fh:
        print(tree.to_str(compact=True, props=['name', 'orf']), file=fh)
    with open(output_prefix + '.node_status.nwk', 'wt') as fh:
        print(tree.write(props=['orf']), file=fh)
    with open(output_prefix + '.named_global_plot.txt', 'wt') as fh:
        print(global_tree.to_str(compact=True, props=['name']), file=fh)
    with open(output_prefix + '.result.txt', 'wt') as fh:
        print(f'origin\t{orf_origin}\tlocal_age\t{local_age}\tglobal_age\t{global_age}', file=fh)
    return orf_origin, local_age, global_age


if __name__ == '__main__':
    infer_orf_age()

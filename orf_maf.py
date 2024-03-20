from ete4 import PhyloTree
import click


def lastn(aln, n):
    i = 0
    while n > 0 and i < len(aln):
        i += 1
        if aln[i] != '-':
            n -= 1
    return i


@click.command(context_settings=dict(
    help_option_names=['-h', '--help'], show_default=True))
@click.argument('bls_path', type=click.File('rt'))
@click.argument('tree_full', type=click.File('rt'))
@click.argument('aln_all', type=click.STRING)
@click.argument('ref_sp', type=click.STRING)
@click.argument('output', type=click.File('wt'))
def orf_maf(bls_path, tree_full, aln_all, ref_sp, output):
    """
    prepare MAF for PhyloCSF calculation (excluding stop codon)

    \b
    bls_path: output of orf_bls.py
    tree_full: the path to global tree with the local_tree as a subtree
    aln_all: alignment of all sequences extracted from whole genome MAF (cactus)
    ref_sp: reference species name in tree/fasta
    output: output file name (could be stdout "-")
    """
    # read the full tree with alignment
    tfull = PhyloTree(tree_full, parser=1, alignment=aln_all, alg_format='fasta')
    # number internal nodes using the same manner as in orf_bls.py
    n = 1
    for node in tfull.traverse('preorder'):
        if not node.is_leaf:
            node.name = f'N{n}'
            n += 1

    # extract origination node
    bls_result = bls_path.readline().split('\t')

    # get leaves under origination node and put ref_sp at the beginning
    origin_age_global = bls_result[5]
    origin_tfull_leaves = [n.name for n in tfull[origin_age_global]]
    if origin_tfull_leaves[0] != ref_sp:
        origin_tfull_leaves.remove(ref_sp)
        origin_tfull_leaves.insert(0, ref_sp)
    naive_age_global = bls_result[17]
    naive_tfull_leaves = [n.name for n in tfull[naive_age_global]]
    if naive_tfull_leaves[0] != ref_sp:
        naive_tfull_leaves.remove(ref_sp)
        naive_tfull_leaves.insert(0, ref_sp)

    # make maf
    # Note that:
    #    - gaps present in all sequences do not influence phylocsf++ results
    #    - size field in s lines do not influence phylocsf++ results
    #    - sequences consisting of gaps only do not influence phylocsf++ results
    ## get start (in reverse order) of start codons
    sv_lastn = lastn(tfull[ref_sp].get_prop('sequence'), 3)
    # initial maf header
    print('track name=ncorf\n##maf version=1\n\na score=100.0', end='\n', file=output)
    # maf block for sequences below origination node determined from ASR
    idlen_max_origin = max([len(n + '.origin') for n in origin_tfull_leaves])
    for n in origin_tfull_leaves:
        name = (n + '.origin').ljust(idlen_max_origin)
        seq = tfull[n].get_prop('sequence')
        if seq:
            seq = seq[:-sv_lastn]
            sline = ' '.join(['s', name, '0', str(len(seq)), '+', '1000000', seq])
            print(sline, end='\n', file=output)
    # maf block for sequences below origination node determined with the naive method
    print('\na score=100.0', end='\n', file=output)
    idlen_max_naive = max([len(n + '.naive') for n in naive_tfull_leaves])
    for n in naive_tfull_leaves:
        name = (n + '.naive').ljust(idlen_max_naive)
        seq = tfull[n].get_prop('sequence')
        if seq:
            seq = seq[:-sv_lastn]
            sline = ' '.join(['s', name, '0', str(len(seq)), '+', '1000000', seq])
            print(sline, end='\n', file=output)
    return


if __name__ == '__main__':
    orf_maf()

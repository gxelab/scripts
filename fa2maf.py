"""
definition of MAF (https://genome.ucsc.edu/FAQ/FAQformat.html#format5):
The first line of a custom MAF track must be a "track" line that contains 
a name=value pair specifying the track name.

The first line of a .maf file begins with ##maf. This word is followed by
white-space-separated variable=value pairs. There should be no white space 
surrounding the "=".
 version - Required. Currently set to one.
 scoring - Optional. A name for the scoring scheme used for the alignments. 

The "s" lines together with the "a" lines define a multiple alignment. The 
first "s" line must be the reference genome. The "s" lines have the following 
fields which are defined by position.
 src -- The name of one of the source sequences for the alignment. For 
         sequences that are resident in a browser assembly, the form 
         'database.chromosome' allows automatic creation of links to other 
         assemblies. Non-browser sequences are typically reference by the species 
         name alone.
 start -- The start of the aligning region in the source sequence. This is a 
         zero-based number. If the strand field is "-" then this is the start 
         relative to the reverse-complemented source sequence (see Coordinate 
         Transforms).
 size -- The size of the aligning region in the source sequence. This number 
         is equal to the number of non-dash characters in the alignment text 
         field below.
 strand -- Either "+" or "-". If "-", then the alignment is to the reverse-
         complemented source.
 srcSize -- The size of the entire source sequence, not just the parts involved
         in the alignment.
 text -- The nucleotides (or amino acids) in the alignment and any insertions 
         (dashes) as well.
"""

from Bio import SeqIO
import click


@click.command(context_settings=dict(
    help_option_names=['-h', '--help'], show_default=True))
@click.argument('fasta', type=click.File('rt'))
@click.argument('maf', type=click.File('wt'))
def fa2maf(fasta, maf):
    """
    convert fasta alignment to a mock MAF(Multiple Alignment Format)
    
    \b
    fasta: input fasta alignment (can be stdin "-")
    maf  : output MAF (can be stdout "-")
    
    note: MAF is for whole genome alignment and contains more information 
    than fasta and it is impossible to convert fasta to a true MAF block.
    This script only generate a pseudo MAF that can be used to calculate
    PhyloCSF score with `phylocsf++` software. The `src` (chromosome),
    `strand` and `srcSize` columns are fill with arbitrary values as place
    holder (0, +, 1000000, respectively). Note that gaps present in all
    sequences, sequences consisting of gaps only, and size field in the s 
    lines do not influence phylocsf++ results.

    \b
    workflow:
    ```
    python fa2maf.py test.aln.fa test.aln.maf
    phylocsf++ score-msa 100vertebrates test.aln.maf
    cat test.aln.maf.scores
    ```
    """
    records = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
    
    print('track name=ncorf\n##maf version=1\n\na score=100.0', end='\n', file=maf)
    
    idlen_max = max([len(sp + '.aln') for sp in records])
    
    for sp in records:
        name = (sp + '.aln').ljust(idlen_max)
        seq = str(records[sp].seq)
        sline = ' '.join(['s', name, '0', str(len(seq)), '+', '1000000', seq])
        print(sline, end='\n', file = maf)
    return


if __name__ == '__main__':
    fa2maf()

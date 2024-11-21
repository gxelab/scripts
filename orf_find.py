# notes:
#  For an ORF with alternative start codons, existing tools like the 
#  start-of-the-art `orfipy` only output the ORF with the most upstream
#  start codon. This script outputs all possible ORFs starting with a
#  given list of start codons and ending with a given list of stop codons,
#  meanning all alternative start codons will be printed with the stop codon
#  as individual ORFs.

from concurrent.futures import ThreadPoolExecutor
import click
from Bio import SeqIO


@click.command(context_settings=dict(help_option_names=['-h', '--help'], show_default=True))
@click.argument('input_file', type=click.Path(exists=True))
@click.option('-b', '--start_codons', default='ATG', help='start codons')
@click.option('-e', '--stop_codons', default='TAG,TAA,TGA', help='stop codons')
@click.option('-d', '--direction', default='both', help='strand direction', type=click.Choice(['fwd', 'rev', 'both']))
@click.option('-t', '--threads', default=1, help='number of threads')
def find_orfs(input_file, start_codons, stop_codons, direction, threads):
    """
    a python script to find all possible open reading frames starting with  a given
    list of start codons and ending with a given list of stop codons from a fasta
    file with multiple sequence entries. 
    output is printed to stdout in tsv format with the following columns:
    seq_id, strand, orf_start, orf_end, orf_len, start_codon, stop_codon, orf_seq
    """
    start_codons = start_codons.split(',')
    stop_codons = stop_codons.split(',')
    
    def process_record(record):
        results = []
        if direction == 'both':
            strands = ['+', '-'] 
        else:
            strands = ['+' if direction == 'fwd' else '-']
        for s in strands:
            for start_codon in start_codons:
                for orf in find_orf(record.seq, start_codon, stop_codons, s):
                    results.append('\t'.join([
                        record.id, s, str(orf[0]), str(orf[1]), 
                        str(orf[1] - orf[0]), start_codon, str(orf[2][-3:]), str(orf[2])]))
        return results

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(process_record, record) for record in SeqIO.parse(input_file, 'fasta')]
        for future in futures:
            for result in future.result():
                print(result)


def find_orf(seq, start_codon, stop_codons, strand):
    start_codon = start_codon.upper()
    stop_codons = [i.upper() for i in stop_codons]
    if strand == '-':
        seq = seq.reverse_complement()
    for i in range(0, len(seq) - 2):
        if seq[i : i + 3] == start_codon:
            for j in range(i + 3, len(seq) - 2, 3):
                if seq[j:j+3] in stop_codons:
                    if strand == '-':
                        yield (len(seq) - (j + 3), len(seq) - i, seq[i : j + 3])
                    else:
                        yield (i, j + 3, seq[i : j + 3])
                    break


if __name__ == '__main__':
    find_orfs()

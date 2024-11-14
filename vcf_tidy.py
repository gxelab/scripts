import sys
import click
from pysam import VariantFile, VariantHeader


@click.command(context_settings=dict(help_option_names=['-h', '--help'], show_default=True))
@click.argument('path_vcf', type=click.STRING)
@click.option('-c', '--chroms', type=click.STRING, default='', help='comma-separated ordered list of chromosomes')
def tidy_main(path_vcf, chroms=''):
    """
    tidy VCF files generated from Drosophila Genome Nexus

    \b
    PATH_VCF: path to the input VCF file (indexed)
    \b
    example:
    ```
    python vcf_tidy.py test.vcf.gz -c chr2L,chr2R
    ```
    """
    vcf_in = VariantFile(path_vcf)
    vcf_contigs = [chrom for chrom in vcf_in.header.contigs]
    # check chroms list
    chroms = chroms.split(',') if chroms != '' else []
    chroms = [chrom for chrom in chroms if chrom in vcf_contigs]
    if chroms == []:
        print('...No valid chrs selected. Using all chrs in the header of the input vcf!', file=sys.stderr)
        chroms = vcf_contigs
    print(f'...chroms selected:{str(chroms)}', file=sys.stderr)
    
    # define a new header
    nh = VariantHeader()
    for chrom in chroms:
        nh.add_line(f'##contig=<ID={chrom}>')
    nh.add_line('##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">')
    nh.add_line('##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">')
    nh.add_line('##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">')
    print(nh, end='')
    
    # parse records from the old vcf
    for chrom in chroms:
        for rec in vcf_in.fetch(chrom):
            if len(rec.alleles) == 1 or rec.alts == ('*',):
                continue
            # substract Simulans from AC & AN
            idx_sim = rec.samples['Simulans']['GT'][0]
            rec.info['AN'] -= 1
            if idx_sim > 0:
                ac = list(rec.info['AC'])
                ac[idx_sim -1] -= 1
                rec.info['AC'] = ac
            # make new record and exclude *
            alleles = [i for i in rec.alleles if i != '*']
            new_rec = nh.new_record(contig=rec.chrom, start=rec.start, stop=rec.stop, alleles=alleles, id=rec.id)
            try:
                alt_n = rec.alts.index('*')
                new_rec.info['AC'] = [j for i, j in enumerate(rec.info['AC']) if i != alt_n]
                new_rec.info['AN'] = rec.info['AN'] - rec.info['AC'][alt_n]
            except ValueError:
                new_rec.info['AC'] = rec.info['AC']
                new_rec.info['AN'] = rec.info['AN']
            if rec.alleles[idx_sim] != '*':
                new_rec.info['AA'] = rec.alleles[idx_sim]
            else:
                new_rec.info['AA'] = 'N'
            print(new_rec, end='')
    return


if __name__ == '__main__':
    tidy_main()


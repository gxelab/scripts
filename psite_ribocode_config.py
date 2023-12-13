import re
import os
import click


@click.command(context_settings=dict(help_option_names=['-h', '--help'], show_default=True))
@click.argument('bam_path', type=click.STRING)
@click.argument('psite_log', type=click.STRING)
@click.option('-p', '--prefix', default='tmp', type=click.STRING,
              help='output prefix')
def get_config(bam_path, psite_log, prefix = 'tmp'):
    """
    get config with *transcriptome.psite.bam
    """
    out_name = f'{prefix}_pre_config.txt'
    outfile = open(out_name, 'w')

    if __name__:
        header = [
            '# SampleName', 'AlignmentFile', 'Stranded(yes/reverse)',
            'P-siteReadLength', 'P-siteLocations']
        print('\t'.join(header), file=outfile)

    alignment_bam = bam_path
    with open(psite_log, 'rt') as fh:
        for line in fh:
            if line.startswith('# Read lengths used:'):
                qwidth_range = re.sub(r'.*?\[(.*?)\].*', r'\1', line.rstrip())
                break
    qwidth_range = re.sub(' ', '', qwidth_range)
    out = [
        os.path.basename(alignment_bam), alignment_bam,
        'yes', qwidth_range, re.sub(r'\d+', '12', qwidth_range)
    ]
    print('\t'.join(out), file=outfile)
    outfile.close()
    
    return

if __name__ == '__main__':
    get_config()

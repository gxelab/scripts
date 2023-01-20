import requests, sys

def get_homolog(gene_id):
    """get all homolog from Ensembl REST API"""
    
    server = "https://rest.ensembl.org"
    ext = "/homology/id/" + gene_id + "?"

    r = requests.get(server + ext, headers = {'Content-Type': 'application/json'})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    
    decoded = r.json()
    for homolog in decoded['data'][0]['homologies']:
        out = [homolog['type'],
               homolog['taxonomy_level'],
               homolog['target']['species'],
               homolog['target']['taxon_id'],
               homolog['target']['id'],
               homolog['target']['protein_id'],
               homolog['target']['perc_pos'],
               homolog['target']['perc_id'],
               homolog['target']['cigar_line'],
               homolog['target']['align_seq'],
               homolog['source']['id'],
               homolog['source']['protein_id'],
               homolog['source']['perc_pos'],
               homolog['source']['perc_id'],
               homolog['source']['cigar_line'],
               homolog['source']['align_seq']]
        print('\t'.join(str(i) for i in out))
    return


if __name__ == '__main__':
    get_homolog(sys.argv[1])


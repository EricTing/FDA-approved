import urllib2
from BeautifulSoup import BeautifulSoup

import pandas as pd

access_num = 'DB01048'

def uniportUrls(access_num):
    my_url = 'http://www.drugbank.ca/drugs/' + access_num
    content = urllib2.urlopen(my_url)
    soup = BeautifulSoup(content.read())
    links = soup('a')
    uniport_urls = [l['href'] for l in links if 'uniprot.org' in l['href']]

    return uniport_urls

uniport_pdb_dat = "../dat/pdbsws_chain.txt"

dset = pd.read_csv(uniport_pdb_dat, sep='\s+', index_col=False, header=None)
pdbs = dset[0]
codes = []
for pdb_code in pdbs:
    codes.append(pdb_code)

codes = set(codes)
for code in codes:
    print code

def readMapping(dat):
    """
    mapping from uniprot access number to pdb id and chain
    """
    entries = [l.rstrip().split() for l in file(dat)]
    acs = set([entry[-1] for entry in entries])

    dset = pd.read_csv(dat, sep='\s+', index_col=False, header=None)
    dset.columns = ['pdb', 'chain', 'ac']

    mapping = {}
    cnt = 0
    for ac in acs:
        if ac != '?':
            my_rows = dset[dset['ac'] == ac]
            pdbs = [' '.join(entry) for entry in my_rows.values[:, 0:2]]
            mapping[ac] = pdbs
            cnt += 1
            print cnt

    return mapping

print readMapping(uniport_pdb_dat)

import urllib2
import cPickle

import pandas as pd

from BeautifulSoup import BeautifulSoup

# access_num = 'DB01048'

def uniprtUrls(access_num):
    my_url = 'http://www.drugbank.ca/drugs/' + access_num
    content = urllib2.urlopen(my_url)
    soup = BeautifulSoup(content.read())
    links = soup('a')
    uniprt_urls = [l['href'] for l in links if 'uniprot.org' in l['href']]

    return uniprt_urls

ifn = "../dat/approved_heavy_atom_num.txt"
dset = pd.read_csv(ifn, sep=' ', index_col=False)

id_urls = {}
for my_id in dset['ID']:
    urls = uniprtUrls(my_id)
    print urls
    id_urls[my_id] = urls

ofn = "../dat/drugid_uniprturls.dat"
with open(ofn, 'wb') as f:
    cPickle.dump(id_urls, f)

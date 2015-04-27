import cPickle
import urllib2

from BeautifulSoup import BeautifulSoup

content = urllib2.urlopen('http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/drugport/GetPage.pl?template=drugindex_pdb.html')

soup = BeautifulSoup(content.read())

links = soup('a')

routes = set([l['href'] for l in links if 'drug_id' in l['href']])

urls = ['http://www.ebi.ac.uk' + route for route in set(routes)]

# print urls

def PdbIds(url):
    print url
    content = urllib2.urlopen(url)
    soup = BeautifulSoup(content.read())
    links = soup('a')
    routes = [l['href'] for l in links
              if '/pdbsum/' in l['href']
              and len(l['href']) < 13]
    pdb_ids = [route.split('/')[-1] for route in routes]
    return pdb_ids


def run():
    drugid_pdb = {}
    for url in urls:
        drug_id = url.split('=')[-1]
        drugid_pdb[drug_id] = PdbIds(url)

    ofn = "../dat/drugid_pdb_embl.dat"
    with open(ofn, 'w') as f:
        cPickle.dump(drugid_pdb, f)

def count():
    '''
    count the number of protein pdbs
    '''
    ifn = "../dat/drugid_pdb_embl.dat"
    with open(ifn, 'r') as f:
        dat = cPickle.load(f)
        cnt = 0
        for key, val in dat.iteritems():
            tmp = []
            for l in val:
                tmp += l
            pdbs = [_.split(' ')[0] for _ in tmp]
            cnt += len(set(pdbs))
        print cnt

# run()

count()

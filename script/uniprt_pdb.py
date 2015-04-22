import cPickle
import pandas as pd


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
    print "total to process", len(acs)
    for ac in acs:
        if ac != '?':
            my_rows = dset[dset['ac'] == ac]
            pdbs = [' '.join(entry) for entry in my_rows.values[:, 0:2]]
            mapping[ac] = pdbs
            cnt += 1
            print "processed", cnt

    return mapping

dset = pd.read_csv(uniprt_pdb_dat, sep='\s+', index_col=False, header=None)
pdbs = dset[0]
codes = []
for pdb_code in pdbs:
    codes.append(pdb_code)

codes = list(set(codes))
# print ','.join(codes)
# for code in codes:
#     print code

uniprt_pdb_dat = "../dat/pdbsws_chain.txt"
mapping = readMapping(uniprt_pdb_dat)
output = open('../dat/uniprt_to_pdb.dat', 'wb')
cPickle.dump(mapping, output)
output.close()

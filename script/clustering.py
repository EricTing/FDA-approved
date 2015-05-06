#!/usr/bin/env python

import os
import cPickle
import pandas as pd

df = pd.DataFrame

CLST_RSLT = "../dat/approved.txt.clusters"
DSET_PATH = "../dat/dat.pkl"

if __name__ == "__main__:":
    dset = pd.read_csv(CLST_RSLT, sep="\t", header=None, names=['DRUGBANK_ID', 'ClusterId'])
    TANI_DIR = "/work/jaydy/working/tanimoto"

    similar_ligs = []
    for drug_id in dset['DRUGBANK_ID']:
        lig_tani_dat = os.path.join(TANI_DIR, drug_id + '.dat')
        if not os.path.exists(lig_tani_dat):
            similar_ligs.append('')
        else:
            with open(lig_tani_dat, 'r') as f:
                dat = cPickle.load(f)
                my_ligs = []
                for lig, tani in dat.iteritems():
                    if float(tani) > 0.9:
                        my_ligs.append(lig)
                similar_ligs.append(' '.join(my_ligs))

    dset['SimilarLigs'] = similar_ligs

    with open(DSET_PATH, 'w') as f:
        cPickle.dump(dset, f)

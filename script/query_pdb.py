#!/usr/bin/env python

'''
given a drugbank id, search the corresponding proteins in the pdb bank
'''

import cPickle
import pandas as pd

uniprt_pdb_dat = "../dat/uniprt_to_pdb.dat"
id_uniprt_ifn = "../dat/drugid_uniprturls.dat"

to_url = {}
with open(id_uniprt_ifn, 'r') as f:
    to_url = cPickle.load(f)

to_pdb = {}
with open(uniprt_pdb_dat, 'r') as f:
    to_pdb = cPickle.load(f)

def query(drug_id):
    urls = to_url[drug_id]
    uniprt_ids = [url.split('/')[-1] for url in urls]
    pdbs = []
    for uniprt_id in uniprt_ids:
        if uniprt_id in to_pdb:
            pdbs.append(to_pdb[uniprt_id])
        else:
            print drug_id, uniprt_id # print the missing ones

    return pdbs


approved_list = "../dat/approved_heavy_atom_num.txt"
dset = pd.read_csv(approved_list, sep=' ', index_col=False)

# drug_id = 'DB00131'
pdbs = {}
for drug_id in dset['ID']:
    pdbs[drug_id] = query(drug_id)

ofn = "../dat/drugid_pdb.dat"
with open(ofn, 'wb') as f:
    cPickle.dump(pdbs, f)

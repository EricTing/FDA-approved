import cPickle
import os
import glob
import pybel
import shutil

from cluster_seq import WORK_DIR, gatherLigIds
from clustering import DSET_PATH

# from cluster_seq import loadSubset, subset_ifn, drug2similarlig_ofn

DAT_DIR = "/work/jaydy/dat/fda_pdb_mb"
INSPECT_DIR = "/work/jaydy/working/Inspect/"
CLUSTER_DIR = "/work/jaydy/working/cluster"


def getPrtMedoids(cluster_result_paths):
    idx2pdb = {}
    for path in cluster_result_paths:
        idx = os.path.basename(path).split('.')[0]
        idx = int(idx)
        my_prts = set()
        lines = [l.rstrip() for l in file(path)]
        for line in lines:
            tokens = line.split()
            if tokens[-1] == '*':
                pdb_id = tokens[-2].split('.')[1][1:]
                my_prts.add(pdb_id)

        idx2pdb[idx] = [_ for _ in my_prts]
    return idx2pdb


with open(DSET_PATH, 'r') as f:
    dset = cPickle.load(f)

cluster_result_paths = glob.glob(CLUSTER_DIR + "/*clustered.clstr")
idx2pdb = getPrtMedoids(cluster_result_paths)

idx2ligs = {}
for idx in set(dset['ClusterId']):
    similar_lig_ids = gatherLigIds(dset, idx)
    prt_ids = [name.split('.')[0] for name in similar_lig_ids]
    if idx in idx2pdb and len(idx2pdb[idx]) > 0:
        intersect_prts = set(idx2pdb[idx]) & set(prt_ids)
        if len(intersect_prts) > 0:
            prt_bnd_ligs = set()
            for lig_id in similar_lig_ids:
                prt_id = lig_id.split('.')[0]
                if prt_id in intersect_prts:
                    prt_bnd_ligs.add(lig_id)
            idx2ligs[idx] = prt_bnd_ligs


prt_bnd_ligs = []
for _, row in dset.iterrows():
    idx = row['ClusterId']
    similar_ligs = set(row['SimilarLigs'].split())
    to_append = ''
    if idx in idx2ligs:
        my_prt_bnd_ligs = similar_ligs & idx2ligs[idx]
        if len(my_prt_bnd_ligs) > 0:
            to_append = ' '.join([_ for _ in my_prt_bnd_ligs])

    prt_bnd_ligs.append(to_append)

dset['ProteinBoundLig'] = prt_bnd_ligs
with open(DSET_PATH, 'w') as f:
    cPickle.dump(dset, f)

# for _, row in dset.iterrows():
#     if row['ProteinBoundLig'] != '':
#         ligs = row['ProteinBoundLig']
#         for lig in ligs.split():
#             print row['DRUGBANK_ID'], lig

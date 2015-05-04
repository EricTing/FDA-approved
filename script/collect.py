import cPickle
import os
import glob
import pybel
import shutil

from cluster_seq import WORK_DIR
from cluster_seq import loadSubset, subset_ifn, drug2similarlig_ofn

DAT_DIR = "/work/jaydy/dat/fda_pdb_mb"
INSPECT_DIR = "/work/jaydy/working/Inspect/"

CLUSTER_DIR = "/work/jaydy/working/cluster"
subsets = loadSubset(subset_ifn)

with open(drug2similarlig_ofn, 'r') as f:
    drug2similarlig = cPickle.load(f)

cluster_result_paths = glob.glob(CLUSTER_DIR + "/*clustered.clstr")


def getPrtMedoids(cluster_result_paths):
    idx2pdb = {}
    for path in cluster_result_paths:
        idx = os.path.basename(path).split('.')[0]
        idx = int(idx)
        if idx in subsets:
            idx2pdb[idx] = set()
            lines = [l.rstrip() for l in file(path)]
            for line in lines:
                tokens = line.split()
                if tokens[-1] == '*':
                    pdb_id = tokens[-2].split('.')[1][1:]
                    idx2pdb[idx].add(pdb_id)

    return idx2pdb


def getMedoidsBoundLigs(idx2pdb):
    drug2lig = {}
    for idx, pdbs in idx2pdb.iteritems():
        drug_ids = subsets[idx]
        for drug_id in drug_ids:
            similar_ligs = drug2similarlig[drug_id]
            if len(similar_ligs) > 0:
                myligs = []
                for lig in similar_ligs:
                    pdb_id = lig.split('.')[0]
                    if pdb_id in pdbs:
                        myligs.append(lig)

                if len(myligs) > 0:
                    drug2lig[drug_id] = myligs
    return drug2lig


idx2pdb = getPrtMedoids(cluster_result_paths)
drug2lig = getMedoidsBoundLigs(idx2pdb)

approved_ifn = "../dat/approved.txt"
drugs = [_ for _ in pybel.readfile("sdf", approved_ifn)]

diff_sz_ofn = "../dat/different_sz.txt"
representative_ofn = "../dat/representative_drugs.txt"
diff_of = open(diff_sz_ofn, 'w')
representative_of = open(representative_ofn, 'w')

tot_drugs = 0
tot_ligs = 0
for drug in drugs:
    drugbank_id = drug.data['DRUGBANK_ID']
    drug.removeh()
    if 8 <= len(drug.atoms) <= 44:
        if drugbank_id in drug2lig:
            tot_drugs += 1
            tot_ligs += len(drug2lig[drugbank_id])
            ligs = drug2lig[drugbank_id]
            for lig in ligs:
                my_id = lig.split('.')[0]
                dat_dir = os.path.join(DAT_DIR, my_id[1:3])
                pdb_path = os.path.join(dat_dir, my_id + '.pdb')
                lig_path = os.path.join(dat_dir, lig)
                assert(os.path.exists(pdb_path))
                assert(os.path.exists(lig_path))

                lig_sz = len([_ for _ in file(lig_path)])
                if lig_sz != len(drug.atoms):
                    diff_of.writelines("%s %s %d %d\n" % (drugbank_id, lig,
                                                          len(drug.atoms), lig_sz))

                representative_of.writelines("%s %s\n" % (drugbank_id, lig))

                inspect_dir = os.path.join(INSPECT_DIR, drugbank_id)
                try:
                    os.makedirs(inspect_dir)
                except:
                    pass
                shutil.copy(pdb_path, inspect_dir)
                shutil.copy(lig_path, inspect_dir)

                drug_ofn = os.path.join(inspect_dir, drugbank_id + ".pdb")
                drug.write("pdb", drug_ofn, overwrite=True)

diff_of.close()
representative_of.close()

print tot_drugs
print tot_ligs

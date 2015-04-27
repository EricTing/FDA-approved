import cPickle
import os

from cluster_seq import WORK_DIR

FASTA_DIR = "/work/jaydy/working/fastas"
FIXED_DIR = "/work/jaydy/working/clean_pdb"

def readDat(ifn):
    with open(ifn, 'r') as f:
        return cPickle.load(f)

dat_ifn = "../dat/drugid_checked_pdb.dat"
dat = readDat(dat_ifn)


new_dat = {}
total = 0
empty = 0
for drug_id, ligs in dat.iteritems():
    size = len(ligs)
    if size == 0:
        empty += 1
    elif size == 1:
        pdb_id = ligs[0].split('.')[0]
        pdb_fn = os.path.join(FIXED_DIR, pdb_id[1:3], pdb_id + '_fix.pdb')
        assert(os.path.exists(pdb_fn))
        with open(pdb_fn, 'r') as f:
            lines = f.readlines()
            seq_len = int(lines[-2].split()[4])
            if seq_len < 600:
                new_dat[drug_id] = ligs
                total += size
            else:
                empty += 1
    else:
        cluster_ifn = os.path.join(FASTA_DIR, drug_id + '.clustered')
        assert(os.path.exists(cluster_ifn))

        my_ligs = set()
        with open(cluster_ifn, 'r') as f:
            lines = f.readlines()
            pdb_ids = [l.split()[0].split('/')[-1] for l in lines if '>' in l]
            for pdb_id in pdb_ids:
                for lig in ligs:
                    if pdb_id in lig:
                        my_ligs.add(lig)

        if len(my_ligs) == 0:
            empty += 1
        new_dat[drug_id] = [i for i in my_ligs]
        total += len(my_ligs)

dat_ofn = "../dat/drugid_clustered_ligs.dat"
with open(dat_ofn, 'w') as f:
    cPickle.dump(new_dat, f)

print total, 'ligand conformations'
print empty, 'drugs with not lig structure found'

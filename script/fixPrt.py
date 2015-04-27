import cPickle
import os
import subprocess32

DAT_DIR = "/work/jaydy/dat/fda_pdb_mb"
WORK_DIR = "/work/jaydy/working/clean_pdb"


def readDat(ifn):
    with open(ifn, 'r') as f:
        return cPickle.load(f)

dat_ifn = "../dat/drugid_checked_pdb.dat"
dat = readDat(dat_ifn)

for drug_id, ligs in dat.iteritems():
    pdb_ids = [name.split('.')[0] for name in ligs]
    for pdb_id in pdb_ids:
        my_work_dir = ""
        try:
            work_dir = os.path.join(WORK_DIR, pdb_id[1:3])
            my_work_dir = work_dir
            os.makedirs(work_dir)
        except:
            pass
        pdb_fn = os.path.join(DAT_DIR, pdb_id[1:3], pdb_id + '.pdb')
        if os.path.exists(pdb_fn):
            os.chdir(my_work_dir)
            ctrip_cmd = ['ctrip', '-prm', '2', pdb_fn]
            subprocess32.call(ctrip_cmd)


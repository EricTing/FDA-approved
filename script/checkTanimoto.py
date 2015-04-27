import cPickle
import os
from count_atms import readSections, CompoundLines

TANI_DIR = "/work/jaydy/working/tanimoto"


ifn = "../dat/approved.txt"
sections = readSections(ifn)

new_dat = {}
for sect in sections:
    cmpdl = CompoundLines(sect)
    drug_id = cmpdl.getDrugBankID()
    try:
        my_dat_ifn = os.path.join(TANI_DIR, drug_id + '.dat')
        with open(my_dat_ifn, 'r') as f:
            dat = cPickle.load(f)
            pdbs = []
            for pdb_id, t_val in dat.iteritems():
                if float(t_val) > 0.9:
                    pdbs.append(pdb_id)
                    print drug_id, pdb_id
    except:
        pass
    new_dat[drug_id] = pdbs

dat_ofn = "../dat/drugid_checked_pdb.dat"
with open(dat_ofn, 'w') as f:
    cPickle.dump(new_dat, f)

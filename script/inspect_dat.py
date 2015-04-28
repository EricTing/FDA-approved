import os
import shutil
import subprocess32
import cPickle

import pandas
df = pandas.DataFrame

from count_atms import CompoundLines, readSections

WORK_DIR = "/work/jaydy/working/Inspect"
LIG_DIR = "/work/jaydy/dat/fda_pdb_mb"

try:
    os.makedirs(WORK_DIR)
except:
    pass

dat_ifn = "../dat/drugid_clustered_ligs.dat"
with open(dat_ifn, 'r') as f:
    dat = cPickle.load(f)


ifn = "../dat/approved.txt"
sects = readSections(ifn)
fda_drugs = {}
for sect in sects:
    com_lines = CompoundLines(sect)
    fda_drugs[com_lines.getDrugBankID()] = com_lines

valid_drugs = []
for drug_id, ligs in dat.iteritems():
    if len(ligs) != 0:
        valid_drugs.append(drug_id)

f = open('../dat/different_sz.txt', 'w')
f.writelines("drug_id lig_id drug_sz lig_sz\n")

f2 = open('../dat/drug_lig.txt', 'w')
f2.writelines("drug_id lig_id\n")

for drug_id in valid_drugs:
    ligs = dat[drug_id]
    for lig in ligs:
        pdb_id = lig.split('.')[0]
        lig_fn = os.path.join(LIG_DIR, pdb_id[1:3], lig)
        assert(os.path.exists(lig_fn))
        lig_sz = len([_ for _ in file(lig_fn)])
        drug_sz = len(fda_drugs[drug_id].getHeavyAtomCoordLines())
        f2.writelines("%s %s\n" % (drug_id, lig))
        if (lig_sz != drug_sz):
            f.writelines("%s %s %d %d\n" % (drug_id, lig, drug_sz, lig_sz))

        lines = fda_drugs[drug_id].copyLines()
        work_dir = os.path.join(WORK_DIR, drug_id)
        try:
            os.makedirs(work_dir)
        except:
            pass
        shutil.copy(lig_fn, work_dir)
        ofn = os.path.join(work_dir, drug_id + '.sdf')
        with open(ofn, 'w') as tmp_f:
            tmp_f.writelines(lines)

        drug_pdb = os.path.join(work_dir, drug_id + '.pdb')
        cmd = ['obabel', '-isdf', ofn, '-opdb', '-O' + drug_pdb]
        subprocess32.call(cmd)
        print ' '.join(cmd)

f.close()
f2.close()

import os
import shutil
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

dat_ifn = "../dat/filtered_drugid_clustered_ligs.dat"
with open(dat_ifn, 'r') as f:
    dat = cPickle.load(f)

sml_drug_dset = pandas.read_csv("../dat/representative_drugs.csv")
sml_drug_ids = [_ for _ in sml_drug_dset['id']]

ifn = "../dat/approved.txt"
sects = readSections(ifn)
fda_drugs = {}
for sect in sects:
    com_lines = CompoundLines(sect)
    fda_drugs[com_lines.getDrugBankID()] = com_lines

f = open('../dat/different_sz.txt', 'w')
f.writelines("drug_id lig_id drug_sz lig_sz\n")
f2 = open('../dat/drug_lig.txt', 'w')
f2.writelines("drug_id lig_id\n")

for drug_id, ligs in dat.iteritems():
    if drug_id in sml_drug_ids:
        drug_sz = fda_drugs[drug_id].heavyAtomNum()
        for lig in ligs:
            lig_id = lig.split('.')[0]
            lig_path = os.path.join(LIG_DIR, lig_id[1:3], lig)
            assert(os.path.exists(lig_path))
            lig_sz = len([_ for _ in file(lig_path)])
            f2.writelines("%s %s\n" % (drug_id, lig_id))

            if lig_sz != drug_sz:
                work_dir = os.path.join(WORK_DIR, drug_id)
                print work_dir
                try:
                    os.makedirs(work_dir)
                except:
                    pass

                shutil.copy(lig_path, work_dir)
                prt_path = os.path.join(LIG_DIR, lig_id[1:3], lig_id + '.pdb')
                shutil.copy(prt_path, work_dir)

                drug_ofn = os.path.join(work_dir, drug_id + '.sdf')
                with open(drug_ofn, 'w') as drug_of:
                    drug_of.writelines(fda_drugs[drug_id].copyLines())
                f.writelines("%s %s %d %d\n" % (drug_id, lig_id, drug_sz, lig_sz))

f.close()
f2.close()

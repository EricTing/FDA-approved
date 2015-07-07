import shelve
import os
import pybel
import cPickle
import pandas as pd
import shutil
df = pd.DataFrame

from clustering import DSET_PATH

with open(DSET_PATH, 'r') as f:
    dset = cPickle.load(f)

drugs = [_ for _ in pybel.readfile("sdf", "../dat/approved.txt")]
ids, sizes = [], []
for drug in drugs:
    drug_id = drug.data['DRUGBANK_ID']
    drug.removeh()
    sz = len(drug.atoms)
    sizes.append(sz)
    ids.append(drug_id)
sz_dset = df([ids, sizes]).T
sz_dset.columns = ['DRUGBANK_ID', 'HeavyAtomNum']
dset = dset.merge(sz_dset)


filter_dset = {'DRUGBANK_ID': [],
               'ProteinBoundLig': [],
               'HeavyAtomNum': [],
               'LigSize': [],
               'LigPath': [],}
DAT_DIR = "/work/jaydy/dat/fda_pdb_mb/"
for _, row in dset.iterrows():
    if row['ProteinBoundLig'] != '':
        ligs = row['ProteinBoundLig']
        drug_id = row['DRUGBANK_ID']
        drug_sz = row['HeavyAtomNum']
        for lig in ligs.split():
            prt_id = lig.split('.')[0]
            lig_path = os.path.join(DAT_DIR, prt_id[1:3], lig)
            assert(os.path.exists(lig_path))
            lig_sz = len([_ for _ in file(lig_path)])
            filter_dset['DRUGBANK_ID'].append(drug_id)
            filter_dset['ProteinBoundLig'].append(lig)
            filter_dset['HeavyAtomNum'].append(drug_sz)
            filter_dset['LigSize'].append(lig_sz)
            filter_dset['LigPath'].append(lig_path)

filter_dset = df(filter_dset)   # all datas
sml_dset = filter_dset[(filter_dset['LigSize'] > 7) & (filter_dset['LigSize'] <= 44)] #  within the range
diff_sz_dset = sml_dset[sml_dset['HeavyAtomNum'] != sml_dset['LigSize']] # those ligs with different size compared with drug
same_sz_dset = sml_dset[sml_dset['HeavyAtomNum'] == sml_dset['LigSize']] # those ligs with different size compared with drug

for _, row in diff_sz_dset.iterrows():
    drug_id = row['DRUGBANK_ID']
    my_drug = filter(lambda d: d.data['DRUGBANK_ID'] == drug_id, drugs)[0]
    drug_dir = os.path.join("/work/jaydy/working/Inspect", drug_id)
    try:
        os.makedirs(drug_dir)
    except:
        pass
    drug_ofn = os.path.join(drug_dir, drug_id + '.pdb')
    my_drug.write(format="pdb", filename=drug_ofn, overwrite=True)

    lig_path = row['LigPath']
    prt_path = os.path.join(os.path.dirname(lig_path),
                            os.path.basename(lig_path).split('.')[0] + '.pdb')
    # import ipdb; ipdb.set_trace()
    assert(os.path.exists(prt_path))
    shutil.copy(lig_path, drug_dir)
    shutil.copy(prt_path, drug_dir)

diff_sz_dset.drop(['LigPath'], axis=1).to_csv("../dat/different_sz.csv", index=False, sep="\t")


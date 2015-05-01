import cPickle
import pandas

df = pandas.DataFrame

from count_atms import CompoundLines, readSections

dat_ifn = "../dat/filtered_drugid_clustered_ligs.dat"
with open(dat_ifn, 'r') as f:
    dat = cPickle.load(f)

valid_drugs = []
for drug_id, ligs in dat.iteritems():
    if len(ligs) != 0:
        valid_drugs.append(drug_id)

ifn = "../dat/approved.txt"
sects = readSections(ifn)
dset = []
for sect in sects:
    cmp_lines = CompoundLines(sect)
    my_id = cmp_lines.getDrugBankID()
    if my_id in valid_drugs:
        dset.append((my_id, len(cmp_lines.getHeavyAtomCoordLines())))

dset = df(dset)
dset.columns = ['id','#HeavyAtom']
within_dset = dset[(dset['#HeavyAtom'] >= 8)
                   & (dset['#HeavyAtom'] <= 44)]

within_dset.to_csv("../dat/representative_drugs.csv")

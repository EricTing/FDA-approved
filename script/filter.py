import cPickle


representatives_ifn = "../dat/drug_representatives.txt"
dat_ifn = "../dat/drugid_clustered_ligs.dat"

with open(dat_ifn, 'r') as f:
    dat = cPickle.load(f)

representatives = [l.rstrip() for l in file(representatives_ifn)]

new_dat = {}
total_lig = 0
for drug_id, ligs in dat.iteritems():
    if drug_id in representatives and len(ligs) > 0:
        total_lig += len(ligs)
        new_dat[drug_id] = ligs

print len(new_dat), total_lig

dat_ofn = "../dat/filtered_drugid_clustered_ligs.dat"
with open(dat_ofn, 'w') as f:
    cPickle.dump(new_dat, f)

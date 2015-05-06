import cPickle
import os

from clustering import DSET_PATH

TANI_DIR = "/work/jaydy/working/tanimoto"
CD_HIT = "/project/michal/apps/cd-hit-v4.6-2012-04-25/cd-hit"
WORK_DIR = "/work/jaydy/working/clean_pdb"
FASTA_DIR = "/work/jaydy/working/fastas"

subset_ifn = "../dat/fda_subset.txt"
drug2similarlig_ofn = "../dat/drug2similarlig.dat"


def prtSeqLen(fasta_fn):
    with open(fasta_fn, 'r') as f:
        line = f.readline()
        seq_len = line.split()[-1]
        return int(seq_len)


def gatherLigIds(dset, idx):
    ligs = dset[dset['ClusterId'] == idx]['SimilarLigs'].values
    lig_ids = ' '.join(ligs)
    lig_set = set(lig_ids.split())
    return lig_set


if __name__ == "__main__":
    with open(DSET_PATH, 'r') as f:
        dset = cPickle.load(f)

    indices = set(dset['ClusterId'])
    for idx in indices:
        lig_set = gatherLigIds(dset, idx)
        pdbs = [name.split('.')[0] for name in lig_set]
        fastas = []
        for pdb_id in pdbs:
            fasta_ifn = os.path.join(WORK_DIR, pdb_id[1:3], pdb_id + '_break.fasta')
            if os.path.exists(fasta_ifn):
                seq_len = prtSeqLen(fasta_ifn)
                if seq_len < 600:
                    with open(fasta_ifn, 'r') as f:
                        fastas.append(f.readlines())

        if len(fastas) > 0:
            fastas_ofn = os.path.join(FASTA_DIR, str(idx) + '.fasta')
            with open(fastas_ofn, 'w') as f:
                for line in fastas:
                    f.writelines(line)

            cmd = [CD_HIT, '-i', fastas_ofn, '-o', str(idx) + '.clustered', '-c 0.7']
            print ' '.join(cmd)

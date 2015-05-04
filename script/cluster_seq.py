import cPickle
import os


TANI_DIR = "/work/jaydy/working/tanimoto"
CD_HIT = "/project/michal/apps/cd-hit-v4.6-2012-04-25/cd-hit"
WORK_DIR = "/work/jaydy/working/clean_pdb"
FASTA_DIR = "/work/jaydy/working/fastas"

subset_ifn = "../dat/fda_subset.txt"
drug2similarlig_ofn = "../dat/drug2similarlig.dat"


def loadSubset(ifn):
    sets = {}
    for line in file(subset_ifn):
        drug_id, set_idx = line.split()
        set_idx = int(set_idx)
        if set_idx not in sets:
            sets[set_idx] = set([drug_id])
        else:
            sets[set_idx].add(drug_id)

    return sets



def findSimilarLigIds(drug_id):
    drug2lig_ifn = os.path.join(TANI_DIR, drug_id + '.dat')
    similar_ligs = []
    if os.path.exists(drug2lig_ifn):
        with open(drug2lig_ifn, 'r') as f:
            lig2tani = cPickle.load(f)

        for lig_id, tani in lig2tani.iteritems():
            if float(tani) > 0.9:
                similar_ligs.append(lig_id)
    else:
        print drug2lig_ifn, 'does not exist'

    return similar_ligs

def prtSeqLen(fasta_fn):
    with open(fasta_fn, 'r') as f:
        line = f.readline()
        seq_len = line.split()[-1]
        return int(seq_len)

if __name__ == "__main__":
    sets = loadSubset(subset_ifn)
    drug2similarlig = {}
    total = 0
    for idx, drugs in sets.iteritems():
        myligs = set()
        for drug_id in drugs:
            similar_ligs = findSimilarLigIds(drug_id)
            drug2similarlig[drug_id] = similar_ligs
            myligs |= set(similar_ligs)
            total += len(myligs)

        if len(myligs) > 0:
            pdbs = [name.split('.')[0] for name in myligs]
            fastas = []
            for pdb_id in pdbs:
                fasta_ifn = os.path.join(WORK_DIR, pdb_id[1:3], pdb_id + '_break.fasta')
                # assert(os.path.exists(fasta_ifn))
                if os.path.exists(fasta_ifn):
                    seq_len = prtSeqLen(fasta_ifn)
                    if seq_len < 600:
                        with open(fasta_ifn, 'r') as f:
                            fastas.append(f.readlines())

                else:
                    pass

            fastas_ofn = os.path.join(FASTA_DIR, str(idx) + '.fasta')
            with open(fastas_ofn, 'w') as f:
                for line in fastas:
                    f.writelines(line)

            cmd = [CD_HIT, '-i', fastas_ofn, '-o', str(idx) + '.clustered', '-c 0.7']
            print ' '.join(cmd)

    with open(drug2similarlig_ofn, 'w') as f:
        cPickle.dump(drug2similarlig, f)

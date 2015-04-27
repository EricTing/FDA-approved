import cPickle
import os

CD_HIT = "/project/michal/apps/cd-hit-v4.6-2012-04-25/cd-hit"
WORK_DIR = "/work/jaydy/working/clean_pdb"
FASTA_DIR = "/work/jaydy/working/fastas"

try:
    os.mkdir(FASTA_DIR)
except:
    pass


def readDat(ifn):
    with open(ifn, 'r') as f:
        return cPickle.load(f)

def prtSeqLen(fasta_fn):
    with open(fasta_fn, 'r') as f:
        line = f.readline()
        seq_len = line.split()[-1]
        return int(seq_len)

if __name__ == "__main__":
    dat_ifn = "../dat/drugid_checked_pdb.dat"
    dat = readDat(dat_ifn)

    print "cd", FASTA_DIR
    for drug_id, ligs in dat.iteritems():
        size = len(ligs)
        if size > 1:
            pdbs = [name.split('.')[0] for name in ligs]
            fastas = []
            for pdb_id in pdbs:
                fasta_ifn = os.path.join(WORK_DIR, pdb_id[1:3], pdb_id + '_break.fasta')
                assert(os.path.exists(fasta_ifn))
                seq_len = prtSeqLen(fasta_ifn)
                if seq_len < 600:
                    with open(fasta_ifn, 'r') as f:
                        fastas.append(f.readlines())

            fastas_ofn = os.path.join(FASTA_DIR, drug_id + '.fasta')
            with open(fastas_ofn, 'w') as f:
                for line in fastas:
                    f.writelines(line)

            os.chdir(FASTA_DIR)
            cmd = [CD_HIT, '-i', fastas_ofn, '-o', drug_id + '.clustered', '-c 0.7']
            print ' '.join(cmd)


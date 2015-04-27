import cPickle
import os

from glob import glob

DAT_DIR = "/work/jaydy/dat/fda_pdb_mb"
WORK_DIR = "/work/jaydy/working/clean_pdb"

pdb2pdb = os.path.abspath("pdb2pdb.pl")
pdb2fasta = os.path.abspath("pdb2fasta.pl")
fasta2fasta = os.path.abspath("fasta2fasta.pl")



def readDat(ifn):
    with open(ifn, 'r') as f:
        return cPickle.load(f)

dat_ifn = "../dat/drugid_checked_pdb.dat"
dat = readDat(dat_ifn)


def run():
    dirs = glob(WORK_DIR + "/*")
    for sub_dir in dirs:
        os.chdir(sub_dir)
        fixed_pdbs = glob("./*_fix.pdb")
        for fixed_pdb in fixed_pdbs:
            pdb_id = fixed_pdb.split('_')[0]
            chain_id = pdb_id[-1]
            cmd0 = ['cd', sub_dir, '&&']
            cmd1 = ['perl', pdb2pdb, '1', chain_id, fixed_pdb, '>',
                    pdb_id + '_shifted.pdb']

            cmd2 = ['perl', pdb2fasta, pdb_id + '_shifted.pdb', pdb_id + '.fasta']

            cmd3 = ['perl', fasta2fasta, '80', pdb_id, pdb_id + '.fasta', pdb_id + '_break.fasta']

            print ' '.join(cmd0 + cmd1)
            print ' '.join(cmd0 + cmd2)
            print ' '.join(cmd0 + cmd3)
            print "echo %s done" % (sub_dir)

if __name__ == "__main__":
    run()

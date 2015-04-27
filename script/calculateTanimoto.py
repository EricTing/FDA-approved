import glob
import os
import cPickle
import subprocess32
from multiprocessing import Pool


from count_atms import readSections, CompoundLines

PDB_DIR = "/work/jaydy/dat/fda_pdb_mb"
WORK_DIR = "/work/jaydy/working/kcombu_run"
TANI_DIR = "/work/jaydy/working/tanimoto"


dat_ifn = "../dat/drugid_pdb.dat"
dat = {}
with open(dat_ifn, 'r') as f:
    dat = cPickle.load(f)

def calculateTanimota(sect):
    try:
        cmpdl = CompoundLines(sect)
        drug_id = cmpdl.getDrugBankID()
        my_pdbs = dat[drug_id]
        ofn = os.path.join(WORK_DIR, drug_id + "_fda.sdf")
        with open(ofn, 'w') as f:
            f.writelines(sect)

        vals = {}
        for pdb_id in my_pdbs:
            pdb_dir = os.path.join(PDB_DIR, pdb_id[1:3])
            my_ligs = glob.glob(pdb_dir + "/*LG*.pdb")
            for lig in my_ligs:
                cmd = ['pkcombu', '-A', ofn, '-B', lig]
                out = subprocess32.Popen(cmd, stdout=subprocess32.PIPE)
                lig_id = os.path.basename(lig)
                tanimoto = [l for l in out.stdout if "#tanimoto" in l][0].split()
                if len(tanimoto) > 1:
                    vals[lig_id] = tanimoto[-1]

        my_dat_ofn = os.path.join(TANI_DIR, drug_id + '.dat')
        with open(my_dat_ofn, 'w') as f:
            cPickle.dump(vals, f)
        print drug_id, "done"
    except:
        print "fails"
        pass

if __name__ == "__main__":
    ifn = "../dat/approved.txt"
    sections = readSections(ifn)
    p = Pool()
    p.map(calculateTanimota, sections)

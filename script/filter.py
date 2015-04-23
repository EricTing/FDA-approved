import errno
import cPickle
import os
import shutil
import subprocess32

WORK_DIR = "/work/jaydy/working/fda_test/"
SCIRPT_DIR = "/home/jaydy/Workspace/Bitbucket/fda-approved-drugs/script"
DAT_DIR = "/work/jaydy/dat/mmCIF"

cif_to_pdb = os.path.join(SCIRPT_DIR, 'qmol_CIF2PDB.pl')
dat_ifn = os.path.join(SCIRPT_DIR, '../dat/drugid_pdb.dat')

def queryMyPdbs(drug_id):
    pdbs = []
    with open(dat_ifn, 'r') as f:
        dat = cPickle.load(f)
        pdbs = dat[drug_id]
        tmp = []
        for l in pdbs:
            tmp += l            # flatten
        pdbs = tmp
    return pdbs

def count():
    '''
    count the number of pdbs
    '''
    with open(dat_ifn, 'r') as f:
        dat = cPickle.load(f)
        cnt = 0
        for key, val in dat.iteritems():
            tmp = []
            for l in val:
                tmp += l
            pdbs = [_.split(' ')[0] for _ in tmp]
            cnt += len(set(pdbs))
        print cnt


my_pdbs = queryMyPdbs('DB05109')
# print my_pdbs
pdb = my_pdbs[0]
pdb_id = pdb.split()[0]
print pdb_id

def runCif2Pdb(work_dir, **kargs):
    os.chdir(work_dir)
    cmd = ['perl', cif_to_pdb, kargs['cif_gz'], pdb_id]
    # os.system(cmd)
    # subprocess32.check_call(cmd)
    print cmd
    subprocess32.call(cmd)

def filterPdb(pdb_id):
    work_dir = os.path.join(WORK_DIR, pdb_id[1:3], pdb_id)
    try:
        os.makedirs(work_dir)
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise e
        pass

    dat_dir = os.path.join(DAT_DIR, pdb_id[1:3], pdb_id)
    cif_gz = os.path.join(DAT_DIR, pdb_id[1:3], pdb_id + '.cif.gz')
    shutil.copy(cif_gz, work_dir)

    cif_gz = pdb_id + '.cif.gz'
    runCif2Pdb(work_dir, cif_gz=cif_gz, pdb_id=pdb_id)

filterPdb(pdb_id)


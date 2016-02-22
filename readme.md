- [FDA-approved small drug 3D structure dataset](#sec-1)
  - [1600 FDA-approved small molecule drugs](#sec-1-1)
  - [Molecule sizes](#sec-1-2)
    - [average heavy atom numbers and standard deviation](#sec-1-2-1)
  - [3D structures of the chemicals](#sec-1-3)
    - [Drug Bank id to pdb id](#sec-1-3-1)
    - [Filter the unbound proteins](#sec-1-3-2)
    - [Cluster the bound proteins based on sequence similarity](#sec-1-3-3)

# FDA-approved small drug 3D structure dataset<a id="orgheadline22"></a>

## 1600 FDA-approved small molecule drugs<a id="orgheadline1"></a>

1.  <http://www.drugbank.ca/downloads#structures>
2.  ./dat/approved.txt

## Molecule sizes<a id="orgheadline5"></a>

### average heavy atom numbers and standard deviation<a id="orgheadline4"></a>

1.  script

    ./script/count\_atms.py

2.  result

    ./dat/approved\_heavy\_atom\_num.txt
    ./dat/elements.txt

## DONE 3D structures of the chemicals<a id="orgheadline21"></a>

### DONE Drug Bank id to pdb id<a id="orgheadline11"></a>

1.  pdb &#x2013;> uniprt

    ./dat/pdbsws\_chain.txt

2.  drug bank id &#x2013;> uniprt ids

    1.  ./script/crawler.py
    2.  ./dat/drugid\_uniprturls.dat

3.  uniprt id &#x2013;> pdb id

    1.  ./script/uniprt\_pdb.py
    2.  ./dat/uniprt\_to\_pdb.dat

4.  drug bank id &#x2013;> pdb id

    1.  ./script/query\_pdb.py
    2.  ./dat/drugid\_pdb.dat
    3.  82773 pdbs
    4.  138078 chains

5.  table

    ```python
    import pickle
    import pandas as pd
    
    dat_ifn = "./dat/drugid_pdb.dat"
    
    with open(dat_ifn, 'r') as f:
        dat = pickle.load(f)
    
    dset = pd.DataFrame({key: list(val) for key, val in dat.iteritems()})
    ```

### DONE Filter the unbound proteins<a id="orgheadline15"></a>

1.  DONE Extract the ligand atoms from the CIF files

    1.  ./script/qmol\_CIF2PDB.pl (not compatible with latest perl modules)
    2.  *work/jaydy/dat/fda\_pdb\_mb*

2.  DONE Cluster drug ligands

    1.  cmd
        
        ```sh
        $ cd ~/Workspace/Bitbucket/fda-approved-drugs/dat/
        $ python /home/jaydy/Workspace/Bitbucket/geauxcluster/src/geauxcluster.py -i approved.txt -o approved.txt.clusters -k DRUGBANK_ID
        ```
    2.  script
        ./script/clustering.py

3.  DONE Compare the pdb ligand with drug bank ligand

    1.  calculate tanimoto coefficient using Kcombu
        1.  ./script/calculateTanimota.py
        2.  pkcombu Segmentation fault for
            1.  DB01049
            2.  to check
                pkcombu -A /work/jaydy/working/kcombu\_run/DB01049\_fda.sdf -B /work/jaydy/dat/fda\_pdb\_mb/ib/3ibdA.LG3.pdb
            3.  DB00707
    2.  refuse proteins if their ligand's tanimoto < 0.9
        1.  ./script/calculateTanimoto.py
        2.  TANI\_DIR = "/work/jaydy/working/tanimoto"

### TODO Cluster the bound proteins based on sequence similarity<a id="orgheadline20"></a>

1.  DONE add missing atoms

    ./script/fixPrt.py

2.  DONE clean the proteins sequences

    python ./script/pdb\_to\_fasta.py > run.sh
    sh ./script/run.sh
    
    1.  ctrip, heavy atom model
        1.  /project/michal/apps/jackal\_64bit/bin/ctrip
    2.  ./script/pdb2pdb.pl
        move pdb sequence to start at 1
    3.  ./script/pdb2fasta.pl
        convert pdb to fasta
    4.  ./script/fasta2fasta.pl
        break the lines at 80

3.  DONE cd-hit to cluster the sequences

    1.  [X] do not count the prt with seq length > 600
    2.  cluster
        python ./script/cluster\_seq.py > ./script/cluster\_seq.sh
        cd /work/jaydy/working/cluster
        sh /home/jaydy/Workspace/Bitbucket/fda-approved-drugs/script/cluster\_seq.sh
    3.  collect
        /home/jaydy/Workspace/Bitbucket/fda-approved-drugs/script/collect.py

4.  TODO 

    1.  original data set of 1554 fda-approved drugs
        1.  script
            ./script/count\_atms.py
            ./dat/drug\_size.txt
        2.  result
            
            <table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">
            
            
            <colgroup>
            <col  class="org-left" />
            
            <col  class="org-right" />
            </colgroup>
            <tbody>
            <tr>
            <td class="org-left">count</td>
            <td class="org-right">1554.000000</td>
            </tr>
            
            
            <tr>
            <td class="org-left">mean</td>
            <td class="org-right">26.257400</td>
            </tr>
            
            
            <tr>
            <td class="org-left">std</td>
            <td class="org-right">18.172654</td>
            </tr>
            
            
            <tr>
            <td class="org-left">min</td>
            <td class="org-right">1.000000</td>
            </tr>
            
            
            <tr>
            <td class="org-left">25%</td>
            <td class="org-right">18.000000</td>
            </tr>
            
            
            <tr>
            <td class="org-left">50%</td>
            <td class="org-right">23.000000</td>
            </tr>
            
            
            <tr>
            <td class="org-left">75%</td>
            <td class="org-right">30.000000</td>
            </tr>
            
            
            <tr>
            <td class="org-left">max</td>
            <td class="org-right">419.000000</td>
            </tr>
            </tbody>
            </table>
        3.  one atom drugs
            
            <table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">
            
            
            <colgroup>
            <col  class="org-left" />
            
            <col  class="org-left" />
            </colgroup>
            <thead>
            <tr>
            <th scope="col" class="org-left">DrugBank ID</th>
            <th scope="col" class="org-left">Chem</th>
            </tr>
            </thead>
            
            <tbody>
            <tr>
            <td class="org-left">DB01356</td>
            <td class="org-left">Li</td>
            </tr>
            
            
            <tr>
            <td class="org-left">DB01370</td>
            <td class="org-left">Al</td>
            </tr>
            
            
            <tr>
            <td class="org-left">DB01592</td>
            <td class="org-left">Fe</td>
            </tr>
            
            
            <tr>
            <td class="org-left">DB01593</td>
            <td class="org-left">Zn</td>
            </tr>
            </tbody>
            </table>
        4.  2 &sigma; range
            (0, 62.6)
        5.  1 &sigma; range
            (8, 44)
    2.  after processing
        1.  script
            ./script/drug\_size.py
        2.  result
            
            <table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">
            
            
            <colgroup>
            <col  class="org-left" />
            
            <col  class="org-right" />
            </colgroup>
            <thead>
            <tr>
            <th scope="col" class="org-left">statistics</th>
            <th scope="col" class="org-right">#HeavyAtom</th>
            </tr>
            </thead>
            
            <tbody>
            <tr>
            <td class="org-left">count</td>
            <td class="org-right">197.000000</td>
            </tr>
            
            
            <tr>
            <td class="org-left">mean</td>
            <td class="org-right">23.928934</td>
            </tr>
            
            
            <tr>
            <td class="org-left">std</td>
            <td class="org-right">12.470170</td>
            </tr>
            
            
            <tr>
            <td class="org-left">min</td>
            <td class="org-right">6.000000</td>
            </tr>
            
            
            <tr>
            <td class="org-left">25%</td>
            <td class="org-right">16.000000</td>
            </tr>
            
            
            <tr>
            <td class="org-left">50%</td>
            <td class="org-right">21.000000</td>
            </tr>
            
            
            <tr>
            <td class="org-left">75%</td>
            <td class="org-right">29.000000</td>
            </tr>
            
            
            <tr>
            <td class="org-left">max</td>
            <td class="org-right">93.000000</td>
            </tr>
            </tbody>
            </table>
        3.  within the range of (8, 44)
            1.  dat
                ./dat/representative\_drugs.csv
                
                <table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">
                
                
                <colgroup>
                <col  class="org-left" />
                
                <col  class="org-right" />
                </colgroup>
                <thead>
                <tr>
                <th scope="col" class="org-left">stats</th>
                <th scope="col" class="org-right">#HeavyAtom</th>
                </tr>
                </thead>
                
                <tbody>
                <tr>
                <td class="org-left">count</td>
                <td class="org-right">186.000000</td>
                </tr>
                
                
                <tr>
                <td class="org-left">mean</td>
                <td class="org-right">22.639785</td>
                </tr>
                
                
                <tr>
                <td class="org-left">std</td>
                <td class="org-right">8.456204</td>
                </tr>
                
                
                <tr>
                <td class="org-left">min</td>
                <td class="org-right">8.000000</td>
                </tr>
                
                
                <tr>
                <td class="org-left">25%</td>
                <td class="org-right">16.000000</td>
                </tr>
                
                
                <tr>
                <td class="org-left">50%</td>
                <td class="org-right">21.000000</td>
                </tr>
                
                
                <tr>
                <td class="org-left">75%</td>
                <td class="org-right">29.000000</td>
                </tr>
                
                
                <tr>
                <td class="org-left">max</td>
                <td class="org-right">44.000000</td>
                </tr>
                </tbody>
                </table>
            2.  3D structures
                1.  protein bound ligands
                    
                    ```python
                    import itertools
                    import shutil
                    import os
                    import pandas as pd
                    import cPickle
                    
                    represents_ifn = "./dat/representative_drugs.csv"
                    dat_ifn = "./dat/dat.pkl"
                    
                    with open(dat_ifn, 'r') as f:
                        dset = cPickle.load(f)
                    
                    represents = pd.read_csv(represents_ifn, index_col=0)
                    
                    inter = dset[dset.DRUGBANK_ID.isin(represents.id)]
                    inter.to_csv("./dat/fda_drugs_prt_lig.csv")
                    
                    drugs = {}
                    for index, row in inter.iterrows():
                        drugs[row.DRUGBANK_ID] = row.ProteinBoundLig.split()
                    
                    complexes = set(itertools.chain(*drugs.values()))
                    print "# complexes %d" % (len(complexes))
                    
                    proteins = set([complex_id.split('.')[0][:4] for complex_id in complexes])
                    print "# proteins %d" % (len(proteins))
                    
                    # retrive the 3D structures
                    DAT_DIR = "/work/jaydy/dat/fda_pdb_mb"
                    STRUCTURE_DIR = "/work/jaydy/dat/fda_drug_structures"
                    
                    def getPaths(ligand_id):
                        mid_two = ligand_id[1:3]
                        ligand_path = os.path.join(DAT_DIR, mid_two, ligand_id)
                        prt_id = ligand_id.split('.')[0] + '.pdb'
                        prt_path = os.path.join(DAT_DIR, mid_two, prt_id)
                        return ligand_path, prt_path
                    
                    def copyFiles(drug_id):
                        dat_dir = os.path.join(STRUCTURE_DIR, drug_id)
                        os.makedirs(dat_dir)
                        for ligand_id in drugs[drug_id]:
                            ligand_path, prt_path = getPaths(ligand_id)
                            shutil.copy(ligand_path, dat_dir)
                            shutil.copy(prt_path, dat_dir)
                    
                        print drug_id, "done"
                    
                    
                    for drug_id in drugs:
                        copyFiles(drug_id)
                    ```

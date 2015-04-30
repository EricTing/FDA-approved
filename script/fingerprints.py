#!/usr/bin/env python

import os
import shutil

from multiprocessing import Pool

import count_atms

WORK_DIR = "/work/jaydy/working/fingers"
FAIL_DIR = "/work/jaydy/working/ff_fail"
TODO_DIR = "/work/jaydy/working/obgen_todo"

try:
    os.makedirs(WORK_DIR)
    os.makedirs(FAIL_DIR)
    os.makedirs(TODO_DIR)
except:
    pass

os.chdir(WORK_DIR)

path = "/work/jaydy/working/fda_sdfs"
cmds = []
for dirpath, dirnames, filenames in os.walk(path):
    for fn in filenames:
        if '3d' in fn:
            drug_id = fn.split('.')[0]
            fullpath = os.path.join(dirpath, fn)
            lines = [l for l in file(fullpath)]
            if len(lines) < 1:
                sdf_fn = os.path.join(dirpath, drug_id + '.sdf')
                with open(sdf_fn, 'r') as f:
                    cmpd = count_atms.CompoundLines(f.readlines())
                    heavy_atoms = cmpd.heavyAtoms()
            elif "WARNING: damped steplength" in lines[0]:
                print fullpath, "WARNING"
                correct_lines = [l for l in lines
                                 if "WARNING" not in l]
                with open(fullpath, 'w') as f:
                    f.writelines(correct_lines)

            if len(lines) > 1:
                drug_id = fn.split('.')[0]
                cmd = ' '.join(['MACCSKeysFingerprints.pl --output FP --CompoundID ', drug_id, fullpath])
                cmds.append(cmd)

# p = Pool()
# p.map(os.system, cmds)
# with open('log', 'w') as f:
#     for cmd in cmds:
#         try:
#             os.system(cmd)
#         except:
#             f.writelines(cmd)
#             pass


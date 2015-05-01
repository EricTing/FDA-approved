#!/usr/bin/env python

import os
import shutil

from multiprocessing import Pool

import count_atms

DAT_DIR = "/work/jaydy/working/fda_3d/"


os.chdir(DAT_DIR)

cmds = []
for dir_name in os.listdir(DAT_DIR):
    sdf = os.path.join(DAT_DIR, dir_name, dir_name + '.sdf')
    if os.path.exists(sdf):
        cmd = ' '.join(['MACCSKeysFingerprints.pl --output FP --CompoundID ', dir_name, sdf])
        cmds.append(cmd)

WORK_DIR = "/work/jaydy/working/prints/"
try:
    os.makedirs(WORK_DIR)
except:
    pass

os.chdir(WORK_DIR)
for cmd in cmds:
    os.system(cmd)

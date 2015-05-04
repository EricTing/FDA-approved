#!/usr/bin/env python
import os
from multi_process import multi_process

paths = [l.rstrip() for l in file("../dat/clean_sdf_paths.txt")]
cmds = []
for path in paths:
    name = os.path.basename(path)
    dir_name = os.path.dirname(path)
    drug_id = name.split('_')[0]
    path_3d = os.path.join(dir_name, drug_id + '_3d.sdf')
    assert(os.path.exists(path_3d))
    path_smi = os.path.join(dir_name, drug_id + '_smi.sdf')
    cmd = ['perl',
           '/home/jaydy/Workspace/Bitbucket/fda-approved-drugs/script/esimdock_sdf',
           '-s', path_3d,
           '-i', drug_id,
           '-c', '-f',
           '-o', path_smi]
    cmds.append(cmd)

multi_process(cmds)


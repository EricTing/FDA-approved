#!/usr/bin/env python

from multi_process import multi_process

paths = [l.rstrip() for l in file("../dat/clean_sdf_paths.txt")]
cmds = [['python', "./fda_drug3d.py", path] for path in paths]
multi_process(cmds)


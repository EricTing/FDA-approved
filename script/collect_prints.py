#!/usr/bin/env python

import os
from glob import glob

DAT_DIR = "/work/jaydy/working/prints"


fpfs = glob(DAT_DIR + "/*MACCSKeysFP.fpf")

ofn = "../dat/fda_prints.txt"
of = open(ofn, 'w')

for fpf in fpfs:
    with open(fpf, 'r') as f:
        lines = f.readlines()
        dat_line = [l for l in lines
                    if not l.startswith('#')]
        if len(dat_line) > 0:
            my_line = dat_line[0]
            my_id = my_line.split()[0][:-1]
            my_prints = my_line.split()[1]
            of.writelines(my_id + " " + my_prints + "\n")

of.close()

#!/usr/bin/env python

import pybel

paths = [l.rstrip() for l in file("../dat/smi_done.txt")]

ofn = "../dat/fda_prints.txt"
of = open(ofn, 'w')
for path in paths:
    drug = pybel.readfile("sdf", path).next()
    prints = drug.data['FINGERPRINT']
    to_write = drug.title + " " + prints + "\n"
    print drug.title, "done"
    of.writelines(to_write)

of.close()

#!/usr/bin/env python
import os
import pybel

from clean_fda_drugs import WORK_DIR


def optCoords(sdf_path):
    drug_id = os.path.basename(sdf_path).split('_')[0]
    drug = pybel.readfile("sdf", sdf_path).next()
    drug.addh()
    drug.localopt()
    drug.removeh()
    drug.data.clear()
    drug.title = drug_id
    ofn = os.path.join(WORK_DIR, drug_id, drug_id + '_3d.sdf')
    drug.write("sdf", ofn, overwrite=True)
    print drug_id, "done"


def main(sdf_path):
    optCoords(sdf_path)

if __name__ == '__main__':
    import sys
    sdf_path = sys.argv[1]
    main(sdf_path)

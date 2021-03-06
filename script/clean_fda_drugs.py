import os
import pybel

from count_atms import readSections, CompoundLines

WORK_DIR = "/work/jaydy/working/fda_3d"


def cleanSdf(WORK_DIR, call=False):
    fda_drugs_ifn = "../dat/approved.txt"
    sections = readSections(fda_drugs_ifn)
    clean_sdf_ofns = []
    for section in sections:
        cmpd = CompoundLines(section)
        cmpd_id = cmpd.getDrugBankID()
        work_dir = os.path.join(WORK_DIR, cmpd_id)
        try:
            os.makedirs(work_dir)
        except:
            pass
        ofn = os.path.join(work_dir, cmpd_id + '_fda.sdf')

        if call:
            with open(ofn, 'w') as f:
                f.writelines(cmpd.copyLines())

        clean_ofn = os.path.join(work_dir, cmpd_id + '_clean.sdf')
        cmd = ['perl', './remove_elements.pl', ofn, clean_ofn]
        if call:
            os.system(' '.join(cmd))
        clean_sdf_ofns.append(clean_ofn)
    return clean_sdf_ofns


if __name__ == "__main__":

    try:
        os.makedirs(WORK_DIR)
    except:
        pass

    clean_sdf_fns = cleanSdf(WORK_DIR, call=True)

    with open("../dat/clean_sdf_paths.txt", 'w') as f:
        for fn in clean_sdf_fns:
            f.writelines(fn + "\n")

    # drugs = []
    # for sdf in clean_sdf_fns:
    #     drug_id = os.path.basename(sdf).split('_')[0]
    #     drug = pybel.readfile("sdf", sdf).next()
    #     drug.title = drug_id
    #     drugs.append(drug)

    # # remove hydrogen atoms
    # for drug in drugs:
    #     drug.removeh()

    # sml_drugs = [drug for drug in drugs
    #             if 8 <= len(drug.atoms) <= 44]

    # for drug in drugs:
    #     optCoords(drug, WORK_DIR)

#!/usr/bin/env python

def readSections(ifn):
    sections = []
    with open(ifn, 'r') as f:
        section = []
        for line in f:
            section.append(line)
            if "$$$$" in line:
                sections.append(section)
                section = []

    return sections


class CompoundLines:
    def __init__(self, section_lines):
        self.lines = [l for l in section_lines]

    def copyLines(self):
        return [_ for _ in self.lines]

    def getID(self):
        return self.lines[1].split()[1]

    def getDrugBankID(self):
        for idx, val in enumerate(self.lines):
            if '<DRUGBANK_ID>' in val:
                return self.lines[idx + 1].rstrip()

    def getAtomNum(self):
        first_token = self.lines[3].split()[0]
        if len(first_token) == 5:
            return int(first_token[:2])
        elif len(first_token) == 6:
            return int(first_token[:3])
        else:
            return int(self.lines[3].split()[0])

    def getAtomCoordLines(self):
        lines = [l for l in self.lines[4:(4 + self.getAtomNum())]]
        return lines

    def getHeavyAtomCoordLines(self):
        return [l for l in self.getAtomCoordLines()
                if l.split()[3] != 'H']

    def heavyAtomNum(self):
        return len(self.getHeavyAtomCoordLines())

    def heavyAtoms(self):
        return [l.split()[3] for l in self.getHeavyAtomCoordLines()]

if __name__ == "__main__":
    ifn = "../dat/approved.txt"
    sects = readSections(ifn)

    elements = set()
    for sect in sects:
        cmp_lines = CompoundLines(sect)
        heavy_atoms = cmp_lines.heavyAtoms()
        elements |= set(heavy_atoms)

    for ele in elements:
        print ele

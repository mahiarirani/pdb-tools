from Bio.PDB import *

parser = PDBParser()
structure = parser.get_structure("2RSC", "2RSC.pdb")
residues = []
for model in structure:
    for chain in model:
        for residue in chain:
            residues.append(residue)


from Bio.PDB import *

def probe(n=100):
    points = []
    phi = math.pi * (3. - math.sqrt(5.))

    for i in range(n):
        y = 1 - (i / float(n - 1)) * 2
        radius = math.sqrt(1 - y * y)
        theta = phi * i
        x = math.cos(theta) * radius
        z = math.sin(theta) * radius

        points.append((x, y, z))

    return points

parser = PDBParser()
structure = parser.get_structure("2RSC", "2RSC.pdb")
residues = []
for model in structure:
    for chain in model:
        for residue in chain:
            residues.append(residue)


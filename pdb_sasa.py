import numpy as np
from Bio.PDB import *
from scipy.spatial import KDTree


#TODO
# - include atoms size
# - output by type
# - find right neighbours
# - input pdb file
# - visualize atom
# - classify the program
# - optimize
# + change probe
# - input probe radius
# - add timer
# - convert percentage to area measurement unit

def probe(n, coords, radius):
    indices = np.arange(0, n, dtype=float) + 0.5

    phi = np.arccos(1 - 2 * indices / n)
    theta = np.pi * (1 + 5 ** 0.5) * indices

    points = np.dstack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)])
    return (points[0] + coords) * radius


def get_coordinates(atom):
    return [atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2]]


coords = []
radius = 1.4
probe_points = 10
parser = PDBParser()
structure = parser.get_structure("2RSC", "2RSC.pdb")
atoms = np.array(Selection.unfold_entities(structure, "A"))
for atom in atoms:
    coords.append(get_coordinates(atom))
    atom.accessibility = 0
    atom.probe = probe(probe_points, atom.get_coord(), radius)

tree = KDTree(coords)
distances, ndx = tree.query(coords, k=10)

for key in range(atoms.size):
    for idx in range(1, 9):
        if distances[key][idx] > radius * 2:
            break
        for point in atoms[ndx[key][0]].probe:
            d = ((point - atoms[ndx[key][idx]].coord) ** 2).sum()
            if d < radius ** 2:
                atoms[ndx[key][0]].accessibility += 1
                break
    print('Atom #%s SASA is %s' % (ndx[key][0], atoms[ndx[key][0]].accessibility / 100))

total = 0
for atom in atoms:
    total += atom.accessibility
print('Total SASA of %s is %s percent' % (structure.get_id(), round(total / atoms.size * 100, 2)))

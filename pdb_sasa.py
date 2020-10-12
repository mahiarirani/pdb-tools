import csv
import numpy as np
from Bio.PDB import *
from sklearn.neighbors import KDTree


#TODO
# + include atoms size
# - output by type
# + find right neighbours
# - input pdb file
# - visualize atom
# - classify the program
# + optimize
# + change probe
# - input arguments
# - add timer
# - convert percentage to area measurement unit

def load_atom_radii(file='atom_radii.csv'):
    atom_radii_dict = {}
    with open(file, 'r') as data:
        for line in csv.DictReader(data):
            atom_radii_dict[line['atom']] = int(line['radii']) / 100
    return atom_radii_dict


def probe(n, coords, radius):
    indices = np.arange(0, n, dtype=float) + 0.5

    phi = np.arccos(1 - 2 * indices / n)
    theta = np.pi * (1 + 5 ** 0.5) * indices

    points = np.dstack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)])
    return (points[0] + coords) * radius


def get_coordinates(atom):
    return [atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2]]


atom_radii = load_atom_radii()
coords = []
probe.radius = 1.4
probe.points = 10
parser = PDBParser()
structure = parser.get_structure("2RSC", "2RSC.pdb")
atoms = np.array(Selection.unfold_entities(structure, "A"))
for atom in atoms:
    coords.append(get_coordinates(atom))
    atom.accessibility = 0
    atom.radius = atom_radii[atom.element]
    atom.probe = probe(probe.points, atom.get_coord(), probe.radius + atom.radius)

coords = np.array(coords)
tree = KDTree(coords)

for index, atom in enumerate(atoms):
    atom.accessibility = sum(
        sum([i != 0 for i in
             [tree.query_radius(atom.probe, atom.radius + probe.radius, count_only=True)]])) / probe.points
    print('Atom #%s [%s] SASA is %s' % (index + 1, atom.element, atom.accessibility))

total = sum(atom.accessibility for atom in atoms)
print('Total SASA of %s is %s percent' % (structure.get_id(), round(total / atoms.size * 100, 2)))
import numpy
from Bio.PDB import *
from scipy.spatial import KDTree


# TODO
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

def probe(n=10):
    indices = numpy.arange(0, n, dtype=float) + 0.5

    phi = numpy.arccos(1 - 2 * indices / n)
    theta = numpy.pi * (1 + 5 ** 0.5) * indices

    points = numpy.dstack([numpy.cos(theta) * numpy.sin(phi), numpy.sin(theta) * numpy.sin(phi), numpy.cos(phi)])
    return points[0]


def get_coordinates(atom):
    return [atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2]]


coords = []
radius = 1.4
parser = PDBParser()
structure = parser.get_structure("2RSC", "2RSC.pdb")
atoms = numpy.array(Selection.unfold_entities(structure, "A"))
for atom in atoms:
    coords.append(get_coordinates(atom))
    atom.accessibility = 0
    atom.probe = probe()
    atom.probe += atom.get_coord()
    atom.probe *= radius

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

total = 0
for atom in atoms:
    total += atom.accessibility
print(total / atoms.size)

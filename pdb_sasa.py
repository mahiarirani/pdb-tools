import math, numpy
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

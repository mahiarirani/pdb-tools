import numpy as np
from sklearn.neighbors import KDTree
from Timer import Timer
from copy import deepcopy

class Probe:
    timer = Timer()

    def __init__(self, atoms, points, radius):
        self.points = points
        self.probe = self.create_probe(self.points)
        self.radius = radius
        self.atoms = atoms

        self.timer.start()
        self.attach_probe()
        self.timer.stop()

        self.timer.start()
        self.coords = self.get_coordinates()
        self.timer.stop()

        self.timer.start()
        self.tree = self.create_tree(self.coords)
        self.timer.stop()

    def get_coordinates(self):
        print('----------\nFetching Probe Coordinates')
        coords = []
        for item in self.atoms:
            for probe in item.probe.coords:
                coords.append(probe)
        print('Probe Coordinates Fetched Successfully\n----------')
        return np.array(coords)

    @staticmethod
    def create_probe(n):
        indices = np.arange(0, n, dtype=float) + 0.5

        phi = np.arccos(1 - 2 * indices / n)
        theta = np.pi * (1 + 5 ** 0.5) * indices

        points = np.dstack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)])
        return points[0]

    def attach_probe(self):
        probe = self.probe
        buried_list = np.stack([[False]] * self.points)
        print('----------\nBegin Probe Attachment')
        for index, atom in enumerate(self.atoms):
            atom.probe = ProbeItem(probe * (self.radius + atom.radius) + atom.get_coord(), atom, buried_list)
            print('Creating Atom #%s [%s] Probe' % (index + 1, atom.element), end='\r')
        print('Probe Attached Successfully\n----------')

    @staticmethod
    def create_tree(item):
        print('----------\nCreating KDTree')
        tree = KDTree(item)
        print('KDTree Created Successfully\n----------')
        return tree

    def get_points_in_atom_probe(self, atoms):
        radius = [(self.radius + atom.radius) - 0.001 for atom in atoms]
        atoms = [[a.get_coord()[0], a.get_coord()[1], a.get_coord()[2]] for a in atoms]
        return self.tree.query_radius(atoms, radius)


class ProbeItem:
    def __init__(self, coords, atom, buried_list):
        self.coords = coords
        self.buried = deepcopy(buried_list)
        self.atom = atom


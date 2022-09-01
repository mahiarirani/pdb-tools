import numpy as np
from sklearn.neighbors import KDTree
from collections import OrderedDict
from copy import deepcopy


class Probe:

    def __init__(self, atoms, points, radius):
        self.points = points
        self.probe = self.create_probe(self.points)
        self.radius = radius
        self.atoms = atoms
        self.attach_probe()
        self.coords = self.get_coordinates()
        self.tree = self.create_tree(self.coords)

    def get_coordinates(self):
        # fetching probe coordinates
        coords = []
        for item in self.atoms:
            for probe in item.probe.coords:
                coords.append(probe)
        return np.array(coords)

    @staticmethod
    def create_probe(n):
        # creating probe points
        indices = np.arange(0, n, dtype=float) + 0.5
        phi = np.arccos(1 - 2 * indices / n)
        theta = np.pi * (1 + 5 ** 0.5) * indices
        points = np.dstack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)])
        return points[0]

    def attach_probe(self):
        # attaching probe points to atoms
        buried_list = np.stack([False] * self.points)
        for index, atom in enumerate(self.atoms):
            atom.probe = _ProbeItem(self.probe * (self.radius + atom.radius) + atom.get_coord(), atom, buried_list)

    @staticmethod
    def create_tree(item):
        # create kd-tree
        return KDTree(item)

    def get_neighbor_probe_points(self, atoms):
        atoms_points = self.get_points_in_atom_probe(atoms)
        for atom, points in enumerate(atoms_points):
            atoms[atom].probe.buried[points % self.points] = True

    def get_points_in_atom_probe(self, atoms, include_probe_radius=True):
        # find points inside an atom probe
        if include_probe_radius:
            radius = [(self.radius + atom.radius) - 0.001 for atom in atoms]
        else:
            radius = [atom.radius for atom in atoms]
        atoms = [a.get_coord() for a in atoms]
        return self.tree.query_radius(atoms, radius)

    def get_points_in_atom_probe_by_distance(self, atoms):
        # find points inside an atom probe sorted by their distance
        radius = [(self.radius + atom.radius) - 0.001 for atom in atoms]
        atoms_coords = [a.get_coord() for a in atoms]
        atoms_points = self.tree.query_radius(atoms_coords, radius, return_distance=True, sort_results=True)
        probes_points = {}
        for index, points in enumerate(atoms_points[0]):
            for idx, point in enumerate(points):
                if probes_points.get(point) is None:
                    probes_points[point] = []
                probes_points[point].append({'atom': index, 'distance': atoms_points[1][index][idx]})
        for item in dict(filter(lambda x: len(x[1]) > 1, probes_points.items())):
            probes_points[item].sort(key=lambda x: x['distance'])
        return atoms_points  # TODO: check return value

    def get_atoms_in_atom_probe(self, atoms, residue=None, atoms_points=None):
        # find atoms inside an atom probe
        if atoms_points is None:
            atoms_points = self.get_points_in_atom_probe(atoms)
        for index, points in enumerate(atoms_points):
            neighbors_points = list(np.array(points) // self.points)
            if residue is None:
                neighbors = [{'atom': self.atoms[atom], 'share': neighbors_points.count(atom)}
                             for atom in list(OrderedDict.fromkeys(neighbors_points))]
            else:
                neighbors = [{'atom': self.atoms[atom], 'share': neighbors_points.count(atom)}
                             for atom in list(OrderedDict.fromkeys(neighbors_points))
                             if self.atoms[atom].get_parent() != residue]
            atoms[index].neighbors = neighbors
        return atoms

    def get_atoms_points_from_neighbors_atoms_probe(self, atoms, residue):
        # find atoms points with neighbors data
        probes_points = [[[] for _ in range(self.points)] for _ in atoms]
        atom_points = [self.get_points_in_atom_probe([a['atom'] for a in atom.neighbors]) for atom in atoms]
        for index, self_atom in enumerate(atom_points):
            for idx, neighbor_atom_points in enumerate(self_atom):
                for i, point in enumerate(neighbor_atom_points):
                    if self.atoms[point // self.points].get_parent() == residue:
                        atom_id = [atom.serial_number for atom in atoms].index(point // self.points + 1)
                        probes_points[atom_id][point % self.points].append(atoms[index].neighbors[idx]['atom'])
        probes_points = [list(f) for f in list(filter(lambda x: len(x) != 0, item) for item in probes_points)]
        for index, self_atom in enumerate(probes_points):
            atoms[index].neighbors = {}
            ratio = (4 * np.pi * (atoms[index].radius + self.radius) ** 2) / self.points
            for atom_neighbors in self_atom:
                for atom in OrderedDict.fromkeys(atom_neighbors):
                    if atoms[index].neighbors.get(atom) is None:
                        atoms[index].neighbors[atom] = 0
                    atoms[index].neighbors[atom] += atom_neighbors.count(atom) / len(atom_neighbors) * ratio
            atoms[index].neighbors = [{'atom': atom, 'share': share} for atom, share in atoms[index].neighbors.items()]
        return atoms


class _ProbeItem:

    def __init__(self, coords, atom, buried_list):
        self.coords = coords
        self.buried = deepcopy(buried_list)
        self.atom = atom

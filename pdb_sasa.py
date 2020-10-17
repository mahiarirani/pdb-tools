import tkinter as tk
from tkinter import filedialog
import os
import csv
import numpy as np
from Bio.PDB import *
from sklearn.neighbors import KDTree


#TODO
# + include atoms size
# + find right neighbours
# + optimize
# + change probe
# + input pdb file
# + classify the program
# + input arguments
# + output by type
# + convert percentage to area measurement unit
# - add timer
# - visualize atom
# - calculate polar/apolar surface
# - concurrent processing

class PDB:
    def __init__(self, probe_points=100, probe_radius=1.4):
        self.probe_points = probe_points
        self.probe_radius = probe_radius

        self.atom_radii = self.load_atom_radii()

        self.structure = self.load_pdb()
        self.atoms = self.get_atoms()
        self.attach_probe()
        self.coords = self.get_all_coordinates()
        self.atom_tree, self.probe_tree = self.create_tree(self.coords)

    @staticmethod
    def load_pdb(file=None):
        if file is None:
            root = tk.Tk()
            root.withdraw()
            file = filedialog.askopenfilename()

        path, file = os.path.split(file)
        file_name, file_extension = os.path.splitext(file)

        parser = PDBParser()
        return parser.get_structure(file_name, path + '/' + file)

    @staticmethod
    def load_atom_radii(file='atom_radii.csv'):
        atom_radii_dict = {}
        with open(file, 'r') as data:
            for line in csv.DictReader(data):
                atom_radii_dict[line['atom']] = int(line['radii']) / 100
        return atom_radii_dict

    def get_atoms(self):
        return Selection.unfold_entities(self.structure, "A")

    @staticmethod
    def get_coordinates(a):
        return [a.get_coord()[0], a.get_coord()[1], a.get_coord()[2]]

    def get_all_coordinates(self):
        atom_coords = []
        probe_coords = []
        for atom in self.atoms:
            atom_coords.append(self.get_coordinates(atom))
            for probe in atom.probe:
                probe_coords.append(probe)
        return {'atom': np.array(atom_coords), 'probe': np.array(probe_coords)}

    @staticmethod
    def create_tree(items):
        trees = []
        print('----------\nBegin Tree Creation')
        for item in items:
            trees.append(KDTree(items[item]))
        print('Tree Created Successfully\n----------')
        return trees

    @staticmethod
    def create_probe(n, atom_coords, atom_radius):
        indices = np.arange(0, n, dtype=float) + 0.5

        phi = np.arccos(1 - 2 * indices / n)
        theta = np.pi * (1 + 5 ** 0.5) * indices

        points = np.dstack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)])
        return points[0] * atom_radius + atom_coords

    def attach_probe(self):
        print('----------\nBegin Probe Attachment')
        for atom in self.atoms:
            atom.radius = self.atom_radii[atom.element]
            atom.probe = self.create_probe(self.probe_points, atom.get_coord(), self.probe_radius + atom.radius)
        print('Probe Attached Successfully\n----------')

    def sasa(self):
        print('----------\nBegin SASA Calculation')
        for index, atom in enumerate(self.atoms):
            radius = (self.probe_radius + atom.radius)
            points, dist = self.probe_tree.query_radius([self.get_coordinates(atom)], radius, return_distance=True)
            for point in points[0]:
                self.coords['probe'][point] = 0
            print('Atom #%s [%s] SASA is calculating' % (index + 1, atom.element))
        for index, atom in enumerate(self.atoms):
            pp = self.probe_points
            atom_probe_points = self.coords['probe'][index * pp:index * pp + pp]
            atom.accessibility = sum([p != 0 for p in atom_probe_points])[0]
            atom.accessibility /= self.probe_points
            atom.accessibility *= 4 * np.pi * (atom.radius + self.probe_radius) ** 2
        print('SASA Calculated Successfully\n----------')

    @staticmethod
    def calc_sasa(item):
        atoms = np.array(Selection.unfold_entities(item, "A"))
        return round(sum(atom.accessibility for atom in atoms), 2)

    def output(self, method=''):
        print('----------\nResult\n----------')
        for model in self.structure:
            if method.__contains__('m'):
                print('Model #%s SASA is %s Å' % (model.id, self.calc_sasa(model)))
            for chain in model:
                if method.__contains__('c'):
                    print('Chain #%s SASA is %s Å' % (chain.id, self.calc_sasa(chain)))
                for residue in chain:
                    if method.__contains__('r'):
                        print('Residue #%s SASA is %s Å' % (residue.get_resname(), self.calc_sasa(residue)))
                    for atom in residue:
                        if method.__contains__('a'):
                            print('Atom #%s SASA is %s Å' % (atom.get_name(), atom.accessibility))
        print('Total SASA of %s is %s Å' % (self.structure.get_id(), self.calc_sasa(self.structure)))


myPDB = PDB()
myPDB.sasa()
myPDB.output()
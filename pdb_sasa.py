import tkinter as tk
from tkinter import filedialog
import os
import csv
import numpy as np
from Bio.PDB import *
from sklearn.neighbors import KDTree
import unknown_radii


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
# + update atom element radii

class PDB:
    def __init__(self, probe_points=100, probe_radius=1.4):
        self.probe_points = probe_points
        self.probe_radius = probe_radius

        self.atom_radii = self.load_atom_radii()

        self.structure = self.load_pdb()
        self.atoms = self.get_atoms()
        self.attach_probe()
        self.coords = np.array(self.get_all_probe_coordinates())
        self.probe_tree = self.create_tree(self.coords)

    @staticmethod
    def load_pdb(file=None):
        print('----------\nLoading PDB File')
        if file is None:
            root = tk.Tk()
            root.withdraw()
            file = filedialog.askopenfilename()

        path, file = os.path.split(file)
        file_name, file_extension = os.path.splitext(file)

        parser = PDBParser()
        structure = parser.get_structure(file_name, path + '/' + file)
        print('PDB File Loaded Successfully\n----------')
        return structure

    @staticmethod
    def load_atom_radii(file='vdw_radii.csv'):
        print('----------\nLoading Atom Radii')
        atom_radii_dict = {}
        with open(file, 'r') as data:
            for line in csv.reader(data):
                if line[0] == 'RESIDUE':
                    residue = line[2]
                    atom_radii_dict[residue] = {}
                elif line[0] == 'ATOM':
                    atom = line[1]
                    atom_radii_dict[residue][atom] = {}
                    atom_radii_dict[residue][atom]['radii'] = float(line[2])
                    atom_radii_dict[residue][atom]['polar'] = bool(line[3])
        print('Atom Radii Loaded Successfully\n----------')
        return atom_radii_dict

    def get_atoms(self):
        return Selection.unfold_entities(self.structure, "A")

    @staticmethod
    def get_coordinates(a):
        return [a.get_coord()[0], a.get_coord()[1], a.get_coord()[2]]

    def get_all_probe_coordinates(self):
        coords = []
        for item in self.atoms:
            for probe in item.probe:
                coords.append(probe)
        return coords

    @staticmethod
    def create_tree(item):
        print('----------\nCreating KDTree')
        tree = KDTree(item)
        print('KDTree Created Successfully\n----------')
        return tree

    @staticmethod
    def create_probe(n):
        indices = np.arange(0, n, dtype=float) + 0.5

        phi = np.arccos(1 - 2 * indices / n)
        theta = np.pi * (1 + 5 ** 0.5) * indices

        points = np.dstack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)])
        return points[0]

    def attach_probe(self):
        probe = self.create_probe(self.probe_points)
        print('----------\nBegin Probe Attachment')
        for index, atom in enumerate(self.atoms):
            res = atom.get_parent().get_resname()
            try:
                atom.radius = self.atom_radii[res][atom.name]['radii']
                atom.polar = self.atom_radii[res][atom.name]['polar']
            except KeyError:
                atom = unknown_radii.get_data(atom)
            atom.probe = probe * (self.probe_radius + atom.radius) + atom.get_coord()
            print('Creating Atom #%s [%s] Probe' % (index + 1, atom.element), end='\r')
        print('Probe Attached Successfully\n----------')

    def sasa(self):
        print('----------\nFinding Near Probe Points')
        for index, atom in enumerate(self.atoms):
            radius = (self.probe_radius + atom.radius)
            points = self.probe_tree.query_radius([self.get_coordinates(atom)], radius)[0]
            self.coords[points] = 0
            print('Searching Atom #%s points' % (index + 1), end='\r')
        print('Probe Points Found Successfully\n----------')
        print('----------\nBegin SASA Calculation')
        for index, atom in enumerate(self.atoms):
            pp = self.probe_points
            atom_probe_points = self.coords[index * pp:index * pp + pp]
            atom.accessibility = sum([p != 0 for p in atom_probe_points])[0]
            atom.accessibility /= pp
            atom.accessibility *= 4 * np.pi * (atom.radius + self.probe_radius) ** 2
            print('Atom #%s [%s] SASA is %s Å' % (index + 1, atom.element, atom.accessibility), end='\r')
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
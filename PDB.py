import tkinter as tk
from tkinter import filedialog
import os
import csv
import unknown_radii
from Bio.PDB import *
from Timer import Timer


class PDB:
    timer = Timer()

    def __init__(self, address=None, atom_radii_file=None, residue_classifications_file=None):
        self.timer.start()
        self.structure = self.load(address)
        self.remove_other_residues()
        self.atoms = self.get_atoms()
        self.timer.stop()

        self.timer.start()
        self.atom_radii = self.load_atom_radii(atom_radii_file)
        self.attach_atom_radii()
        self.timer.stop()

        self.timer.start()
        self.residue_classifications = self.load_residue_classifications(residue_classifications_file)
        self.attach_residue_classification()
        self.timer.stop()

    @staticmethod
    def load(file=None):
        print('----------\nLoading PDB File', end='\r')
        if file is None:
            root = tk.Tk()
            root.withdraw()
            file = filedialog.askopenfilename()

        path, file = os.path.split(file)
        file_name, file_extension = os.path.splitext(file)

        parser = PDBParser()
        parser.QUIET = True
        structure = parser.get_structure(file_name, path + '/' + file)
        print('PDB File Loaded Successfully')
        return structure

    @staticmethod
    def load_atom_radii(file=None):
        print('----------\nLoading Atom Radii', end='\r')
        if file is None:
            file = 'vdw_radii.csv'
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
                    atom_radii_dict[residue][atom]['polar'] = bool(int(line[3]))
        print('Atom Radii Loaded Successfully')
        return atom_radii_dict

    def attach_atom_radii(self):
        for atom in self.atoms:
            res = atom.get_parent().get_resname()
            try:
                atom.radius = self.atom_radii[res][atom.name]['radii']
                atom.polar = self.atom_radii[res][atom.name]['polar']
            except KeyError:
                atom.radius, atom.polar = unknown_radii.get_data(atom.element, atom.name)
        return self.atoms

    @staticmethod
    def load_residue_classifications(file=None):
        print('----------\nLoading Residue Classes', end='\r')
        if file is None:
            file = 'residue_classifications.csv'
        residue_classes_dict = {}
        with open(file, 'r') as data:
            reader = csv.reader(data)
            next(reader)
            for line in reader:
                residue_classes_dict[line[0]] = {}
                residue_classes_dict[line[0]]['polar'] = bool(line[1])
                residue_classes_dict[line[0]]['charge'] = line[1]
        print('Residue Classes Loaded Successfully')
        return residue_classes_dict

    def attach_residue_classification(self):
        residues = self.structure.get_residues()
        for residue in residues:
            try:
                residue.polar = self.residue_classifications[residue.get_resname()]['polar']
                residue.charge = self.residue_classifications[residue.get_resname()]['charge']
            except KeyError:
                residue.polar = False
                residue.charge = 0

    def get_atoms(self, item=None):
        if item is None:
            item = self.structure
        return Selection.unfold_entities(item, "A")

    @staticmethod
    def get_coordinates(a):
        return [a.get_coord()[0], a.get_coord()[1], a.get_coord()[2]]

    def get_item(self, model, chain=None, residue=None):
        try:
            if chain is None:
                item = self.structure[model]
            elif residue is None:
                item = self.structure[model][chain]
            else:
                item = self.structure[model][chain][(' ', residue, ' ')]
        except KeyError:
            item = None
        return item

    def remove_other_residues(self):
        delete_residues = []
        residues = [r for r in self.structure.get_residues()]
        for residue in residues:
            if residue.get_id()[0] != ' ':
                delete_residues.append(residue)
        [delete_residue.get_parent().detach_child(delete_residue.id) for delete_residue in delete_residues]

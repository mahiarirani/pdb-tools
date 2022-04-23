import os
import csv
import unknown_radii
from Bio.PDB import *
from Timer import Timer


class PDB:
    timer = Timer()

    def __init__(self, address, atom_radii_file, residue_classifications_file, residue_rsa_file):
        self.timer.start()
        self.structure = self.load(address)
        self.remove_other_residues()
        self.timer.stop()

        self.timer.start()
        self.atom_radii = self.load_atom_radii(atom_radii_file)
        self.attach_atom_radii()
        self.timer.stop()

        self.timer.start()
        self.residue_classifications = self.load_residue_classifications(residue_classifications_file)
        self.attach_residue_classification()
        self.timer.stop()

        self.timer.start()
        self.residue_RSA = self.load_residue_rsa(residue_rsa_file)
        self.attach_residue_rsa()
        self.timer.stop()

    @staticmethod
    def load(address):
        print('----------\nLoading PDB File', end='\r')
        
        path, file = os.path.split(address)
        file_name, file_extension = os.path.splitext(file)

        parser = PDBParser()
        parser.QUIET = True
        structure = parser.get_structure(file_name, address)
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
        for atom in self.get_atoms():
            res = atom.get_parent().get_resname()
            try:
                atom.radius = self.atom_radii[res][atom.name]['radii']
                atom.polar = self.atom_radii[res][atom.name]['polar']
            except KeyError:
                atom.radius, atom.polar = unknown_radii.get_data(atom.element, atom.name)

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
                residue_classes_dict[line[0]]['polar'] = bool(int(line[1]))
                residue_classes_dict[line[0]]['charge'] = int(line[2])
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

    @staticmethod
    def load_residue_rsa(file=None):
        print('----------\nLoading Residue RSA', end='\r')
        if file is None:
            file = 'relative_solvent_accessibility.csv'
        residue_rsa_dict = {}
        with open(file, 'r') as data:
            reader = csv.reader(data)
            next(reader)
            for line in reader:
                residue_rsa_dict[line[0]] = int(line[1])
        print('Residue rsa Loaded Successfully')
        return residue_rsa_dict

    def attach_residue_rsa(self):
        residues = self.structure.get_residues()
        for residue in residues:
            try:
                residue.RSA = self.residue_RSA[residue.get_resname()]
            except KeyError:
                residue.RSA = 1

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

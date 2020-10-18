import tkinter as tk
from tkinter import filedialog
import os
from Bio.PDB import *


class PDB:
    def __init__(self, address=None):
        self.structure = self.load(address)
        self.atoms = self.get_atoms()

    @staticmethod
    def load(file=None):
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

    def get_atoms(self, item=None):
        if item is None:
            item = self.structure
        return Selection.unfold_entities(item, "A")

    @staticmethod
    def get_coordinates(a):
        return [a.get_coord()[0], a.get_coord()[1], a.get_coord()[2]]


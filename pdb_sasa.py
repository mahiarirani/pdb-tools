import numpy as np
from sklearn.neighbors import KDTree
from Probe import Probe
from PDB import PDB


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


class Main:
    def __init__(self, probe_points=100, probe_radius=1.4):
        self.PDB = PDB()
        self.atoms = self.PDB.get_atoms()

        self.probe = Probe(self.atoms, probe_points, probe_radius)
        self.probe.coords = np.array(self.probe.get_coordinates())
        self.probe.tree = self.create_tree(self.probe.coords)

    @staticmethod
    def create_tree(item):
        print('----------\nCreating KDTree')
        tree = KDTree(item)
        print('KDTree Created Successfully\n----------')
        return tree

    def sasa(self):
        print('----------\nFinding Near Probe Points')
        for index, atom in enumerate(self.atoms):
            radius = (self.probe.radius + atom.radius) - 0.001
            points = self.probe.tree.query_radius([self.PDB.get_coordinates(atom)], radius)[0]
            self.probe.coords[points] = 0
            print('Searching Atom #%s points' % (index + 1), end='\r')
        print('Probe Points Found Successfully\n----------')
        print('----------\nBegin SASA Calculation')
        for index, atom in enumerate(self.atoms):
            pp = self.probe.points
            atom_probe_points = self.probe.coords[index * pp:index * pp + pp]
            atom.accessibility = sum([p != 0 for p in atom_probe_points])[0] / pp
            atom.accessibility *= 4 * np.pi * (atom.radius + self.probe.radius) ** 2
            print('Atom #%s [%s] SASA is %s Å' % (index + 1, atom.element, atom.accessibility), end='\r')
        print('SASA Calculated Successfully\n----------')

    def calc_sasa(self, item):
        atoms = np.array(self.PDB.get_atoms(item))
        return round(sum(atom.accessibility for atom in atoms), 2)

    def report(self, method=''):
        print('----------\nResult\n----------')
        for model in self.PDB.structure:
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
        print('Total SASA of %s is %s Å' % (self.PDB.structure.get_id(), self.calc_sasa(self.PDB.structure)))


myPDB = Main()
myPDB.sasa()
myPDB.report()
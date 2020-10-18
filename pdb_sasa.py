import numpy as np
from Probe import Probe
from PDB import PDB
from Timer import Timer


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
# + add timer
# - visualize atom
# + calculate polar/apolar surface
# - concurrent processing
# + update atom element radii


class Main:
    timer = Timer()

    def __init__(self, probe_points=100, probe_radius=1.4):
        self.PDB = PDB()
        self.probe = Probe(self.PDB.atoms, probe_points, probe_radius)

    def sasa(self):
        self.timer.start()
        print('----------\nFinding Near Probe Points')
        for index, atom in enumerate(self.PDB.atoms):
            radius = (self.probe.radius + atom.radius) - 0.001
            points = self.probe.tree.query_radius([self.PDB.get_coordinates(atom)], radius)[0]
            self.probe.coords[points] = 0
            print('Searching Atom #%s points' % (index + 1), end='\r')
        print('Probe Points Found Successfully\n----------')
        self.timer.stop()
        self.timer.start()
        print('----------\nBegin SASA Calculation')
        for index, atom in enumerate(self.PDB.atoms):
            pp = self.probe.points
            atom_probe_points = self.probe.coords[index * pp:index * pp + pp]
            atom.accessibility = sum([p != 0 for p in atom_probe_points])[0] / pp
            atom.accessibility *= 4 * np.pi * (atom.radius + self.probe.radius) ** 2
            print('Atom #%s [%s] SASA is %s Å' % (index + 1, atom.element, atom.accessibility), end='\r')
        print('SASA Calculated Successfully\n----------')
        self.timer.stop()
        self.timer.lapsed()

    def calc_sasa(self, item):
        polar, apolar = 0, 0
        for atom in self.PDB.get_atoms(item):
            if atom.polar:
                polar += atom.accessibility
            else:
                apolar += atom.accessibility
        polar = round(polar, 2)
        apolar = round(apolar, 2)
        total = round(polar + apolar, 2)
        return total, polar, apolar

    def report(self, method=''):
        print('----------\nResult\n----------')
        for model in self.PDB.structure:
            if method.__contains__('m'):
                t, p, a = self.calc_sasa(model)
                print('Model #%s SASA is %s Å [Polar : %s Å / Apolar : %s Å]' % (
                    model.id, t, p, a))
            for chain in model:
                if method.__contains__('c'):
                    t, p, a = self.calc_sasa(chain)
                    print('Chain #%s SASA is %s Å [Polar : %s Å / Apolar : %s Å]' % (
                        chain.id, t, p, a))
                for residue in chain:
                    if method.__contains__('r'):
                        t, p, a = self.calc_sasa(residue)
                        print('Residue #%s SASA is %s Å [Polar : %s Å / Apolar : %s Å]' % (
                            residue.get_resname(), t, p, a))
                    for atom in residue:
                        if method.__contains__('a'):
                            print('Atom #%s SASA is %s Å [Polar : %s]' % (
                                atom.get_name(), atom.accessibility, atom.polar))
        t, p, a = self.calc_sasa(self.PDB.structure)
        print('Total SASA of %s is %s Å [Polar: %s Å / Apolar: %s Å]' % (
            self.PDB.structure.get_id(), t, p, a))


myPDB = Main()
myPDB.sasa()
myPDB.report()
import numpy as np
from Probe import Probe
from PDB import PDB
from Timer import Timer


class FastSASA:
    timer = Timer()

    def __init__(self, probe_points=100, probe_radius=1.4):
        self.PDB = PDB()
        self.probe = Probe(self.PDB.atoms, probe_points, probe_radius)

    def sasa(self, report=''):
        self.get_neighbor_probe_points()
        self.calc_sasa()
        self.report_sasa(report)

    def get_neighbor_probe_points(self):
        self.timer.start()
        print('----------\nFinding Neighbor Probe Points')
        atoms_points = self.probe.get_points_in_atom_probe(self.PDB.atoms)
        for atom, points in enumerate(atoms_points):
            self.probe.atoms[atom].probe.buried[points % self.probe.points] = True
            print('Searching in Atom #%s Radius' % (atom + 1), end='\r')
        print('Neighbor Probe Points Found Successfully\n----------')
        self.timer.stop()

    def calc_sasa(self):
        self.timer.start()
        print('----------\nBegin SASA Calculation')
        for index, atom in enumerate(self.PDB.atoms):
            pp = self.probe.points
            atom_probe_points = self.probe.atoms[index].probe.buried
            atom.accessibility = sum([not p[0] for p in atom_probe_points]) / pp
            atom.accessibility *= 4 * np.pi * (atom.radius + self.probe.radius) ** 2
            print('Atom #%s [%s] SASA is %s Å' % (index + 1, atom.element, atom.accessibility), end='\r')
        print('SASA Calculated Successfully\n----------')
        self.timer.stop()

    def sum_sasa(self, item):
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

    def report_sasa(self, method=''):
        print('----------\nResult\n----------')
        for model in self.PDB.structure:
            if method.__contains__('m'):
                t, p, a = self.sum_sasa(model)
                print('Model #%s SASA is %s Å [Polar : %s Å / Apolar : %s Å]' % (
                    model.id, t, p, a))
            for chain in model:
                if method.__contains__('c'):
                    t, p, a = self.sum_sasa(chain)
                    print('Chain #%s SASA is %s Å [Polar : %s Å / Apolar : %s Å]' % (
                        chain.id, t, p, a))
                for residue in chain:
                    if method.__contains__('r'):
                        t, p, a = self.sum_sasa(residue)
                        print('Residue #%s SASA is %s Å [Polar : %s Å / Apolar : %s Å]' % (
                            residue.get_resname(), t, p, a))
                    for atom in residue:
                        if method.__contains__('a'):
                            print('Atom #%s SASA is %s Å [Polar : %s]' % (
                                atom.get_name(), atom.accessibility, atom.polar))
        t, p, a = self.sum_sasa(self.PDB.structure)
        print('Total SASA of %s is %s Å [Polar: %s Å / Apolar: %s Å]' % (
            self.PDB.structure.get_id(), t, p, a))

    def residue_neighbors(self, model, chain, residue):
        neighbors = self.get_residue_neighbors(model, chain, residue)
        self.report_neighbors(neighbors)

    def get_residue_neighbors(self, model, chain, residue):
        item = self.PDB.get_item(model, chain, residue)
        if item is None:
            return None
        atoms = self.PDB.get_atoms(item)
        atoms = self.get_atoms_neighbors(atoms)
        print('----------\nBegin Residue Neighbor Search')
        item.neighbors = {}
        for atom in atoms:
            for neighbor in atom.neighbors:
                res = neighbor.get_parent()
                chain = res.get_parent()
                model = chain.get_parent()
                if res != item:
                    if item.neighbors.get(model.id) is None:
                        item.neighbors[model.id] = {}
                    if item.neighbors[model.id].get(chain.id) is None:
                        item.neighbors[model.id][chain.id] = []
                    if res.id[1] not in item.neighbors[model.id][chain.id]:
                        item.neighbors[model.id][chain.id].append(res.id[1])
        print('Residue Neighbors Found Successfully\n----------')
        return item.neighbors

    def get_atoms_neighbors(self, atoms):
        print('----------\nBegin Atom Neighbor Search')
        atoms_points = self.probe.get_points_in_atom_probe(atoms)
        for index, atom_points in enumerate(atoms_points):
            atoms[index].neighbors = [self.probe.atoms[atom] for atom in atom_points // self.probe.points]
        print('Atom Neighbors Found Successfully\n----------')
        return atoms

    def report_neighbors(self, neighbors):
        print('----------\nResult\n----------')
        for model in neighbors:
            print('Model #%s : ' % model)
            for chain in neighbors[model]:
                print('Chain %s : [' % chain, end='')
                for residue in neighbors[model][chain]:
                    print('%s #%s ,' % (self.PDB.get_item(model, chain, residue).get_resname(), residue),
                          end='')
                print('\b\b]')


if __name__ == '__main__':
    myPDB = FastSASA()
    myPDB.sasa()
    myPDB.timer.lapsed()
    myPDB.timer.reset()
    myPDB.residue_neighbors(0, 'A', 20)
    myPDB.timer.lapsed()

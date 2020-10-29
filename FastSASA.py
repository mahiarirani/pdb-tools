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
        print('----------\nFinding Neighbor Probe Points', end='\r')
        atoms_points = self.probe.get_points_in_atom_probe(self.PDB.atoms)
        for atom, points in enumerate(atoms_points):
            self.probe.atoms[atom].probe.buried[points % self.probe.points] = True
            print('Searching in Atom #%s Radius' % (atom + 1), end='\r')
        print('Neighbor Probe Points Found Successfully')
        self.timer.stop()

    def calc_sasa(self):
        self.timer.start()
        print('----------\nBegin SASA Calculation', end='\r')
        for index, atom in enumerate(self.PDB.atoms):
            pp = self.probe.points
            atom_probe_points = self.probe.atoms[index].probe.buried
            atom.accessibility = sum([not p[0] for p in atom_probe_points]) / pp
            atom.accessibility *= 4 * np.pi * (atom.radius + self.probe.radius) ** 2
            print('Atom #%s [%s] SASA is %s Å' % (index + 1, atom.element, atom.accessibility), end='\r')
        print('SASA Calculated Successfully')
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
        print('----------\nResult :\n')
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
        print('Total SASA of %s is %s Å [Polar: %s Å / Apolar: %s Å]\n' % (
            self.PDB.structure.get_id(), t, p, a))

    def residue_neighbors(self, model, chain, residue):
        item = self.PDB.get_item(model, chain, residue)
        if item is None:
            return None
        neighbors = self.get_residue_neighbors(item)
        self.report_residue_neighbors(item, neighbors)

    def get_residue_neighbors(self, residue, quiet=False):
        atoms = self.PDB.get_atoms(residue)
        atoms = self.get_atoms_neighbors(atoms, quiet)
        if not quiet:
            self.timer.start()
            print('----------\nBegin Residue Neighbor Search', end='\r')
        residue.neighbors = {}
        for atom in atoms:
            for neighbor_residue in list(set(n.get_parent() for n in atom.neighbors if n.get_parent)):
                chain = neighbor_residue.get_parent()
                model = chain.get_parent()
                if neighbor_residue != residue:
                    if residue.neighbors.get(model.id) is None:
                        residue.neighbors[model.id] = {}
                    if residue.neighbors[model.id].get(chain.id) is None:
                        residue.neighbors[model.id][chain.id] = []
                    residue.neighbors[model.id][chain.id].append(neighbor_residue.id[1])
                    residue.neighbors[model.id][chain.id] = list(set(residue.neighbors[model.id][chain.id]))
        if not quiet:
            print('Residue Neighbors Found Successfully')
            self.timer.stop()
        return residue.neighbors

    def get_atoms_neighbors(self, atoms, quiet=False):
        if not quiet:
            self.timer.start()
            print('----------\nBegin Atom Neighbor Search', end='\r')
        atoms_points = self.probe.get_points_in_atom_probe(atoms)
        for index, atom_points in enumerate(atoms_points):
            atoms[index].neighbors = [self.probe.atoms[atom] for atom in atom_points // self.probe.points]
        if not quiet:
            print('Atom Neighbors Found Successfully')
            self.timer.stop()
        return atoms

    def report_residue_neighbors(self, item, neighbors):
        print('----------\nResult :\n')
        print('Selected Residue is %s #%s\n' % (item.get_resname(), item.id[1]))
        for model in neighbors:
            print('Model #%s : ' % model)
            for chain in neighbors[model]:
                print('Chain %s : [' % chain, end='')
                for residue in neighbors[model][chain]:
                    print('%s #%s ,' % (self.PDB.get_item(model, chain, residue).get_resname(), residue),
                          end='')
                print('\b\b]')
        print('')

    def chain_neighbors(self, model, chain):
        neighbors = self.get_chain_neighbors(model, chain)
        self.report_chain_neighbors(model, chain, neighbors)

    def get_chain_neighbors(self, model, chain):
        self.timer.start()
        print('----------\nSearching Chain Neighbors', end='\r')
        chain = self.PDB.structure[model][chain]
        chain.neighbors = {}
        for residue in chain.get_residues():
            print('Getting Residue #%s Neighbors' % residue.id[1], end='\r')
            neighbors = self.get_residue_neighbors(residue, True)[0]
            if len(neighbors) > 1:
                ch = residue.get_parent()
                for neighbor in [key for key in neighbors.keys() if key not in [ch.id]]:
                    if ch.neighbors.get(neighbor) is None:
                        ch.neighbors[neighbor] = {}
                    if len(neighbors[neighbor]):
                        ch.neighbors[neighbor][neighbors[neighbor][0]] = True
        for neighbor in chain.neighbors:
            chain.neighbors[neighbor] = list(chain.neighbors[neighbor].keys())
        print('Chain Neighbors Found Successfully')
        self.timer.stop()
        return chain.neighbors

    def report_chain_neighbors(self, model, chain, neighbors):
        print('----------\nResult :\n')
        print('Selected Chain is %s \n' % chain)
        for chain in neighbors:
            print('Chain %s :  [' % chain, end='')
            for residue in neighbors[chain]:
                print('%s #%s ,' % (self.PDB.get_item(model, chain, residue).get_resname(), residue),
                      end='')
            print('\b\b]')
        print('')


if __name__ == '__main__':
    myPDB = FastSASA()
    myPDB.timer.lapsed()
    myPDB.timer.reset()
    myPDB.sasa()
    myPDB.timer.lapsed()
    myPDB.timer.reset()
    myPDB.residue_neighbors(0, 'A', 20)
    myPDB.timer.lapsed()
    myPDB.timer.reset()
    myPDB.chain_neighbors(0, 'A')
    myPDB.timer.lapsed()

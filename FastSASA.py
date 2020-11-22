import numpy as np
from collections import OrderedDict
from Probe import Probe
from PDB import PDB
from Timer import Timer


class FastSASA:
    timer = Timer()

    def __init__(self, probe_points=100, probe_radius=1.4):
        self.PDB = PDB()
        self.probe = Probe(self.PDB.get_atoms(), probe_points, probe_radius)

    def sasa(self, report=''):
        self.get_neighbor_probe_points()
        self.calc_sasa()
        self.sum_sasa(self.PDB.structure)
        self.report_sasa(report)

    def get_neighbor_probe_points(self):
        self.timer.start()
        print('----------\nFinding Neighbor Probe Points', end='\r')
        atoms_points = self.probe.get_points_in_atom_probe(self.PDB.get_atoms())
        for atom, points in enumerate(atoms_points):
            self.probe.atoms[atom].probe.buried[points % self.probe.points] = True
            print('Searching in Atom #%s Radius' % (atom + 1), end='\r')
        print('Neighbor Probe Points Found Successfully')
        self.timer.stop()

    def calc_sasa(self):
        self.timer.start()
        print('----------\nBegin SASA Calculation', end='\r')
        for index, atom in enumerate(self.PDB.get_atoms()):
            pp = self.probe.points
            atom_probe_points = self.probe.atoms[index].probe.buried
            atom.accessibility = sum([not p for p in atom_probe_points]) / pp
            atom.size = 4 * np.pi * (atom.radius + self.probe.radius) ** 2
            atom.sasa = atom.accessibility * atom.size
            print('Atom #%s [%s] SASA is %s Å' % (index + 1, atom.element, atom.sasa), end='\r')
        print('SASA Calculated Successfully')
        self.timer.stop()

    def sum_sasa(self, item):
        if hasattr(item, 'sasa'):
            return item.sasa, item.size, item.accessibility
        else:
            total_sasa, total_size = 0, 0
            for i in item.get_list():
                item_sasa, item_size, _ = self.sum_sasa(i)
                total_sasa += item_sasa
                total_size += item_size
            item.sasa = total_sasa
            item.size = item.RSA if hasattr(item, 'RSA') else total_size
            item.accessibility = item.sasa / item.size * 100
            return item.sasa, item.size, item.accessibility

    def report_sasa(self, method=''):
        print('----------\nResult :\n')
        for model in self.PDB.structure:
            if method.__contains__('m'):
                t, _, a = map(round, self.sum_sasa(model), [2] * 3)
                print('Model #%s SASA is %s Å (%s%%)' % (
                    model.id, t, a))
            for chain in model:
                if method.__contains__('c'):
                    t, _, a = map(round, self.sum_sasa(chain), [2] * 3)
                    print('Chain #%s SASA is %s Å (%s%%)' % (
                        chain.id, t, a))
                for residue in chain:
                    if method.__contains__('r'):
                        t, _, a = map(round, self.sum_sasa(residue), [2] * 3)
                        print('Residue %s #%s SASA is %s Å (%s%%)' % (
                            residue.get_resname(), residue.get_id()[1], t, a))
                    for atom in residue:
                        if method.__contains__('a'):
                            print('Atom #%s SASA is %s Å [Polar : %s]' % (
                                atom.get_name(), atom.sasa, atom.polar))
        t, _, a = map(round, self.sum_sasa(self.PDB.structure), [2] * 3)
        print('Total SASA of %s is %s Å (%s%%)\n' % (
            self.PDB.structure.get_id(), t, a))

    def residue_neighbors(self, model, chain, residue):
        item = self.PDB.get_item(model, chain, residue)
        if item is None:
            print('----------\nError Getting Residue Neighbors :\nResidue not Found\n')
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
            neighbor_residues = [n.get_parent() for n in atom.neighbors if n.get_parent() is not residue]
            for neighbor_residue in list(OrderedDict.fromkeys(neighbor_residues)):
                chain = neighbor_residue.get_parent()
                model = chain.get_parent()
                if neighbor_residue != residue:
                    if residue.neighbors.get(model.id) is None:
                        residue.neighbors[model.id] = {}
                    if residue.neighbors[model.id].get(chain.id) is None:
                        residue.neighbors[model.id][chain.id] = []
                    residue.neighbors[model.id][chain.id].append(neighbor_residue.id[1])
        for model in residue.neighbors:
            for chain in residue.neighbors[model]:
                residue.neighbors[model][chain] = list(OrderedDict.fromkeys(residue.neighbors[model][chain]))
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
        item = self.PDB.get_item(model, chain)
        if item is None:
            print('----------\nError Getting Chain Neighbors :\nChain not Found\n')
            return None
        neighbors = self.get_chain_neighbors(item)
        self.report_chain_neighbors(model, chain, neighbors)

    def get_chain_neighbors(self, chain):
        self.timer.start()
        print('----------\nSearching Chain Neighbors', end='\r')
        chain.neighbors = {}
        for residue in chain.get_residues():
            print('Getting Residue #%s Neighbors' % residue.id[1], end='\r')
            neighbors = self.get_residue_neighbors(residue, True)[chain.get_parent().id]
            if neighbors is not None and len(neighbors) > 1:
                ch = residue.get_parent()
                for neighbor in [key for key in neighbors.keys() if key is not ch.id]:
                    if ch.neighbors.get(neighbor) is None:
                        ch.neighbors[neighbor] = {}
                    if len(neighbors[neighbor]):
                        ch.neighbors[neighbor][residue.id[1]] = neighbors[neighbor][0]
        for neighbor in chain.neighbors:
            temp = {val: key for key, val in chain.neighbors[neighbor].items()}
            chain.neighbors[neighbor] = {val: key for key, val in temp.items()}
        print('Chain Neighbors Found Successfully')
        self.timer.stop()
        return chain.neighbors

    def report_chain_neighbors(self, model, selected_chain, neighbors):
        print('----------\nResult :\n')
        print('Selected Chain is %s \n' % selected_chain)
        for chain in neighbors:
            counter = 0 if len(neighbors[chain]) > 3 else 1
            print('Chain %s :  [' % chain, end='')
            for own_residue in neighbors[chain]:
                if counter % 3 == 0:
                    print('', end='\n\t')
                counter += 1
                neighbor = neighbors[chain][own_residue]
                own_residue_name = self.PDB.get_item(model, selected_chain, own_residue).get_resname()
                neighbor_name = self.PDB.get_item(model, chain, neighbor).get_resname()
                print('%s #%s -> %s #%s ,' % (own_residue_name, own_residue, neighbor_name, neighbor), end='')
            print('\b\b]')
        print('')

    def critical_residues(self, threshold, model, chain):
        if not hasattr(self.PDB.structure, 'sasa'):
            self.sasa()
        item = self.PDB.get_item(model, chain)
        if item is None:
            print('----------\nError Getting Chain Neighbors :\nChain not Found\n')
            return None
        critical_residues = self.get_critical_residues(threshold, model, chain, item)
        self.report_critical_residues(critical_residues, chain)

    def get_critical_residues(self, threshold, model, chain, item):
        self.timer.start()
        print('----------\nSearching For Critical Residues', end='\r')
        critical_residues = []
        residues = item.get_residues()
        for residue in residues:
            if residue.accessibility > threshold:
                continue
            print('Checking Residue #%s Neighbors' % residue.id[1], end='\r')
            neighbors = self.get_residue_neighbors(residue, True)
            equal, non = [], []
            for r in neighbors[model][chain]:
                res = self.PDB.get_item(model, chain, r)
                if res.polar == residue.polar:
                    equal.append(r)
                else:
                    non.append(r)
            if residue.polar:
                if len(equal) == 0:
                    critical_residues.append({'residue': residue, 'type': 'HydPhl-HydPhb', 'neighbors': non})
                else:
                    critical_residues.append({'residue': residue, 'type': 'HydPhl-HydPhl', 'neighbors': equal})
            else:
                if len(non) == 0:
                    critical_residues.append({'residue': residue, 'type': 'HydPhb-HydPhb', 'neighbors': equal})
        print('Critical Residues Found Successfully')
        self.timer.stop()
        return critical_residues

    @staticmethod
    def report_critical_residues(items, selected_chain):
        print('----------\nResult :\n')
        print('Selected Chain is %s \n' % selected_chain)
        items_types_list = [i['type'] for i in items]
        print('%s Critical Residues Found : %s HydPhl-HydPhl - %s HydPhl-HydPhb - %s HydPhb-HydPhb\n' % (
             len(items),
             items_types_list.count('HydPhl-HydPhl'),
             items_types_list.count('HydPhl-HydPhb'),
             items_types_list.count('HydPhl-HydPhb')))

        for item in items:
            residue = item['residue']
            print('Residue %s #%s with %s %% RSA is %s with %s charge' % (
                residue.get_resname(),
                residue.get_id()[1],
                round(residue.accessibility, 2),
                'Hydrophilic' if residue.polar else 'Hydrophobic',
                'Positive' if residue.charge == 1 else 'Negative' if residue.charge == -1 else 'Natural'))
            print('There are %s neighbors : %s\n' % (item['type'], item['neighbors']))
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
    myPDB.timer.reset()
    myPDB.critical_residues(5, 0, 'A')
    myPDB.timer.lapsed()

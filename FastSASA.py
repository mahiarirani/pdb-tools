import numpy as np
import argparse
import pathlib
import pickle
import os
from collections import OrderedDict
from Probe import Probe
from PDB import PDB
from Timer import Timer


class FastSASA:
    timer = Timer()

    def __init__(self, pdb_file, probe_points=100, probe_radius=1.4, atom_radii_file=None,
                 residue_classification_file=None, residue_radii_file=None):
        self.PDB = PDB(pdb_file, atom_radii_file, residue_classification_file, residue_radii_file)
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
                            print('Atom #%s SASA is %s Å' % (
                                atom.get_name(), atom.sasa))
        t, _, a = map(round, self.sum_sasa(self.PDB.structure), [2] * 3)
        print('Total SASA of %s is %s Å (%s%%)\n' % (
            self.PDB.structure.get_id(), t, a))

    def residue_neighbors(self, model, chain, residue):
        item = self.PDB.get_item(model, chain, residue)
        if item is None:
            print('----------\nError Getting Residue Neighbors :\nResidue not Found\n')
            return None
        self.get_residue_neighbors(item)
        self.report_residue_neighbors(item)

    def get_residue_neighbors(self, residue, quiet=False):
        if not quiet:
            self.timer.start()
            print('----------\nBegin Residue Neighbor Search', end='\r')
        atoms = self.PDB.get_atoms(residue)
        atoms = self.probe.get_atoms_in_atom_probe(atoms, residue)
        atoms = self.probe.get_atoms_points_from_neighbors_atoms_probe(atoms, residue)
        residue.neighbors = {}
        for atom in atoms:
            neighbor_residues = {}
            for n in atom.neighbors:
                share = n['share']
                if neighbor_residues.get(n['atom'].get_parent()) is not None:
                    share += neighbor_residues[n['atom'].get_parent()]
                neighbor_residues.update({n['atom'].get_parent(): share})
            for neighbor_residue in neighbor_residues:
                chain = neighbor_residue.get_parent()
                model = chain.get_parent()
                if residue.neighbors.get(model) is None:
                    residue.neighbors[model] = {}
                if residue.neighbors[model].get(chain) is None:
                    residue.neighbors[model][chain] = {}
                share = neighbor_residues[neighbor_residue]
                if residue.neighbors[model][chain].get(neighbor_residue) is not None:
                    share += residue.neighbors[model][chain][neighbor_residue]
                residue.neighbors[model][chain][neighbor_residue] = share
        if not quiet:
            print('Residue Neighbors Found Successfully')
            self.timer.stop()
        return residue.neighbors

    def report_residue_neighbors(self, item):
        print('----------\nResult :\n')
        print('Selected Residue is %s #%s\n' % (item.get_resname(), item.id[1]))
        for model in item.neighbors:
            print('Model #%s : ' % model.id)
            for chain in item.neighbors[model]:
                print('Chain %s : ' % chain.id)
                for residue in sorted(item.neighbors[model][chain]):
                    share = round(item.neighbors[model][chain][residue], 2)
                    print('%s #%s (%s Å)' % (residue.get_resname(), residue.get_id()[1], share))
                print('')
        print('')

        atoms = self.PDB.get_atoms(item)
        for atom in atoms:
            print('%s of %s#%s :' % (atom, item.get_resname(), item.get_id()[1]))
            for n in atom.neighbors:
                neighbor_residue = n['atom'].get_parent()
                print('%s Å contact with %s of %s#%s in Chain %s' % (
                    round(n['share'], 2), n['atom'],
                    neighbor_residue.get_resname(),
                    neighbor_residue.get_id()[1],
                    neighbor_residue.get_full_id()[2]))
            print('')
        print('')

    def chain_neighbors(self, model, chain):
        item = self.PDB.get_item(model, chain)
        if item is None:
            print('----------\nError Getting Chain Neighbors :\nChain not Found\n')
            return None
        self.get_chain_neighbors(item)
        self.report_chain_neighbors(item)

    def get_chain_neighbors(self, chain):
        self.timer.start()
        print('----------\nSearching Chain Neighbors', end='\r')
        chain.neighbors = {}
        for residue in chain.get_residues():
            print('Getting Residue #%s Neighbors' % residue.id[1], end='\r')
            neighbors = self.get_residue_neighbors_shallow(residue, True)[chain.get_parent()]
            if neighbors is not None and len(neighbors) > 1:
                ch = residue.get_parent()
                for neighbor in [key for key in neighbors.keys() if key is not ch]:
                    if ch.neighbors.get(neighbor) is None:
                        ch.neighbors[neighbor] = {}
                    if len(neighbors[neighbor]):
                        ch.neighbors[neighbor][residue] = sorted(neighbors[neighbor])
        print('Chain Neighbors Found Successfully')
        self.timer.stop()
        return chain.neighbors

    def get_residue_neighbors_shallow(self, residue, quiet=False):
        if not quiet:
            self.timer.start()
            print('----------\nBegin Residue Neighbor Search', end='\r')
        atoms = self.PDB.get_atoms(residue)
        atoms = self.probe.get_atoms_in_atom_probe(atoms, residue)
        residue.neighbors = {}
        for atom in atoms:
            neighbor_residues = [n['atom'].get_parent() for n in atom.neighbors]
            for neighbor_residue in neighbor_residues:
                chain = neighbor_residue.get_parent()
                model = chain.get_parent()
                if residue.neighbors.get(model) is None:
                    residue.neighbors[model] = {}
                if residue.neighbors[model].get(chain) is None:
                    residue.neighbors[model][chain] = []
                residue.neighbors[model][chain].append(neighbor_residue)
        for model in residue.neighbors:
            for chain in residue.neighbors[model]:
                residue.neighbors[model][chain] = list(OrderedDict.fromkeys(residue.neighbors[model][chain]))
        if not quiet:
            print('Residue Neighbors Found Successfully')
            self.timer.stop()
        return residue.neighbors

    @staticmethod
    def report_chain_neighbors(item):
        print('----------\nResult :\n')
        print('Selected Chain is %s \n' % item.id)
        for chain in item.neighbors:
            print('Chain %s :' % chain.id)
            for own_residue in item.neighbors[chain]:
                residue_neighbors = [n.get_id()[1] for n in item.neighbors[chain][own_residue]]
                print('%s #%s -> %s ,' % (own_residue.get_resname(), own_residue.get_id()[1], residue_neighbors))
        print('')

    def critical_residues(self, threshold, model, chain):
        if not hasattr(self.PDB.structure, 'sasa'):
            self.sasa()
        item = self.PDB.get_item(model, chain)
        if item is None:
            print('----------\nError Getting Chain Neighbors :\nChain not Found\n')
            return None
        critical_residues = self.get_critical_residues(threshold, item)
        self.report_critical_residues(critical_residues, chain)

    def get_critical_residues(self, threshold, item):
        self.timer.start()
        print('----------\nSearching For Critical Residues', end='\r')
        critical_residues = []
        residues = item.get_residues()
        for residue in residues:
            chain = residue.get_parent()
            model = chain.get_parent()
            if residue.get_id() == 66:
                print(residue.accessibility)
            if residue.accessibility > threshold:
                continue
            print('Checking Residue #%s Neighbors' % residue.id[1], end='\r')
            neighbors = self.get_residue_neighbors_shallow(residue, True)
            equal, non = [], []
            for res in neighbors[model][chain]:
                if res.polar == residue.polar:
                    equal.append(res)
                else:
                    non.append(res)
            if residue.polar:
                if len(equal) == 0:
                    critical_residues.append({'residue': residue, 'type': 'HydPhl-HydPhb', 'neighbors': non})
                else:
                    critical_residues.append({'residue': residue, 'type': 'HydPhl-HydPhl', 'neighbors': equal})
            else:
                if len(non) == 0:
                    critical_residues.append({'residue': residue, 'type': 'HydPhb-HydPhb', 'neighbors': equal})
                else:
                    critical_residues.append({'residue': residue, 'type': 'HydPhb-HydPhl', 'neighbors': non})
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
            print('There are %s neighbors : %s\n' % (item['type'], [n.get_id()[1] for n in item['neighbors']]))
        print('')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process PDB files \n'
                                                 'Abbreviations: M - model, C - chain, R - residue, T - threshold')
    parser.add_argument('pdb', help="pdb file location", type=pathlib.Path)

    group_actions = parser.add_argument_group('actions')
    group_actions.add_argument('-s', '--sasa',
                               help="calculate pbs sasa", action='store_true')
    group_actions.add_argument('-sr', '--sasa-report',
                               help="report sasa in detail", metavar='report', default='')
    group_actions.add_argument('-r', '--residue-neighbors',
                               help="find residue neighbors", nargs=3, metavar=('M', 'C', 'R'))
    group_actions.add_argument('-c', '--chain-neighbors',
                               help="find chain neighbors", nargs=2, metavar=('M', 'C'))
    group_actions.add_argument('-cr', '--critical-residues',
                               help="find critical residue", nargs=3, metavar=('T', 'M', 'C'))

    group_options = parser.add_argument_group('options')
    group_options.add_argument('-pp', '--probe_points',
                               help="number of points on atom water probe", metavar='probe_points', default=100)
    group_options.add_argument('-pr', '--probe_radius',
                               help="radius of atom water probe", metavar='probe_radius', default=1.4)

    group_assets = parser.add_argument_group('assets')
    group_assets.add_argument('-far', '--file-atom-radii',
                              help="atom radii File", type=pathlib.Path, metavar='file_location')
    group_assets.add_argument('-frc', '--file-residue-classification',
                              help="residue classification file", type=pathlib.Path, metavar='file_location')
    group_assets.add_argument('-frr', '--file-residue-rsa',
                              help="file residue rsa file", type=pathlib.Path, metavar='file_location')
    args = parser.parse_args()

    path, file = os.path.split(args.pdb)
    if not os.path.isdir(os.getcwd() + '/PDBObject/'):
        os.makedirs(os.getcwd() + '/PDBObject/')
    if os.path.isfile(os.getcwd() + '/PDBObject/' + file + 'o'):
        with open(os.getcwd() + '/PDBObject/' + file + 'o', 'rb') as f:
            myPDB = pickle.load(f)
    else:
        myPDB = FastSASA(args.pdb, args.probe_points, args.probe_radius, args.file_atom_radii,
                         args.file_residue_classification, args.file_residue_rsa)
        with open(os.getcwd() + '/PDBObject/' + myPDB.PDB.structure.get_id() + '.pdbo', 'wb') as f:
            pickle.dump(myPDB, f)

    myPDB.timer.lapsed()

    if (args.sasa):
        myPDB.timer.reset()
        myPDB.sasa(args.sasa_report)
        myPDB.timer.lapsed()

    if (args.residue_neighbors):
        params = args.residue_neighbors
        myPDB.timer.reset()
        myPDB.residue_neighbors(int(params[0]), params[1], int(params[2]))
        myPDB.timer.lapsed()

    if (args.chain_neighbors):
        params = args.chain_neighbors
        myPDB.timer.reset()
        myPDB.chain_neighbors(int(params[0]), params[1])
        myPDB.timer.lapsed()

    if (args.critical_residues):
        params = args.critical_residues
        myPDB.timer.reset()
        myPDB.critical_residues(int(params[0]), int(params[1]), params[2])
        myPDB.timer.lapsed()

import numpy as np
import argparse
import pathlib
import pickle
import os
from PDBObjectFileManager import FileManager
from collections import OrderedDict
from Probe import Probe
from PDB import PDB, Residue, Chain


class PDBTools:

    def __init__(self, pdb_file, probe_points, probe_radius):
        self.pdb = PDB(pdb_file, probe_points, probe_radius)
        self.fm = FileManager(self.pdb)
        self.fm.ready()

    def sasa(self, report=''):
        self.get_neighbor_probe_points()
        self.calc_sasa()
        self.sum_sasa(self.pdb.structure)
        self.fm.add(['sasa'])
        self.report_sasa(report)

    def get_neighbor_probe_points(self):
        print('----------\nFinding Neighbor Probe Points', end='\r')
        atoms_points = self.pdb.probe.get_points_in_atom_probe(self.pdb.get_atoms())
        for atom, points in enumerate(atoms_points):
            self.pdb.probe.atoms[atom].probe.buried[points % self.pdb.probe.points] = True
            print('Searching in Atom #%s Radius' % (atom + 1), end='\r')
        print('Neighbor Probe Points Found Successfully')

    def calc_sasa(self):
        print('----------\nBegin SASA Calculation', end='\r')
        for index, atom in enumerate(self.pdb.get_atoms()):
            pp = self.pdb.probe.points
            atom_probe_points = self.pdb.probe.atoms[index].probe.buried
            atom.accessibility = sum([not p for p in atom_probe_points]) / pp
            atom.size = 4 * np.pi * (atom.radius + self.pdb.probe.radius) ** 2
            atom.sasa = atom.accessibility * atom.size
            print('Atom #%s [%s] SASA is %s Å' % (index + 1, atom.element, atom.sasa), end='\r')
        print('SASA Calculated Successfully')

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
            item.size = item.RSA if hasattr(item, 'RSA') and item.RSA != 1 else total_size
            item.accessibility = item.sasa / item.size * 100
            return item.sasa, item.size, item.accessibility

    def report_sasa(self, method=''):
        print('----------\nResult :\n')
        for model in self.pdb.structure:
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
        t, _, a = map(round, self.sum_sasa(self.pdb.structure), [2] * 3)
        print('Total SASA of %s is %s Å (%s%%)\n' % (
            self.pdb.structure.get_id(), t, a))

    def residue_neighbors(self, residue: str):
        residue = Residue(residue)
        item = self.pdb.get_item(residue.model, residue.chain, residue.number)
        if item is None:
            print('----------\nError Getting Residue Neighbors :\nResidue not Found\n')
            return None
        self.get_residue_neighbors(item)
        self.fm.add(['neighbors'], item)

    def get_residue_neighbors(self, residue):
        atoms = self.pdb.get_atoms(residue)
        atoms = self.pdb.probe.get_atoms_in_atom_probe(atoms, residue)
        atoms = self.pdb.probe.get_atoms_points_from_neighbors_atoms_probe(atoms, residue)
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
                if residue.neighbors.get(model.get_id()) is None:
                    residue.neighbors[model.get_id()] = {}
                if residue.neighbors[model.get_id()].get(chain.get_id()) is None:
                    residue.neighbors[model.get_id()][chain.get_id()] = {}
                share = neighbor_residues[neighbor_residue]
                if residue.neighbors[model.get_id()][chain.get_id()].get(neighbor_residue.get_id()[1]) is not None:
                    share += residue.neighbors[model.get_id()][chain.get_id()][neighbor_residue.get_id()[1]]
                residue.neighbors[model.get_id()][chain.get_id()][neighbor_residue.get_id()[1]] = share
        for atom in atoms:
            neighbors = {}
            for n in atom.neighbors:
                id = n['atom'].get_full_id()
                c, r, a = id[2], id[3][1], id[4][0]
                if neighbors.get(c) is None:
                    neighbors[c] = {}
                if neighbors[c].get(r) is None:
                    neighbors[c][r] = {}
                if neighbors[c][r].get(a) is None:
                    neighbors[c][r][a] = {}
                neighbors[c][r][a] = n['share']
            atom.neighbors = neighbors
        return residue.neighbors

    def report_residue_neighbors(self, item):
        print('----------\nResult :\n')
        print('Selected Residue is %s #%s\n' % (item.get_resname(), item.id[1]))
        for model in item.neighbors:
            print('Model #%s : ' % model)
            for chain in item.neighbors[model]:
                print('Chain %s : ' % chain)
                for residue in sorted(item.neighbors[model][chain]):
                    share = round(item.neighbors[model][chain][residue], 2)
                    print('%s #%s (%s Å)' % (residue, residue, share))
                print('')
        print('')

        atoms = self.pdb.get_atoms(item)
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

    def chain_neighbors(self, chain):
        chain = Chain(chain)
        item = self.pdb.get_item(chain.model, chain.chain)
        if item is None:
            print('----------\nError Getting Chain Neighbors :\nChain not Found\n')
            return None
        self.get_chain_neighbors(item)
        self.fm.add(['neighbors'], item)

    def get_chain_neighbors(self, chain):
        print('----------\nSearching Chain Neighbors', end='\r')
        chain.neighbors = {}
        for residue in chain.get_residues():
            print('Getting Residue #%s Neighbors' % residue.id[1], end='\r')
            neighbors = self.get_residue_neighbors_shallow(residue)[chain.get_parent()]
            if neighbors is not None and len(neighbors) > 1:
                ch = residue.get_parent()
                for neighbor in [key for key in neighbors.keys() if key is not ch]:
                    if ch.neighbors.get(neighbor.get_id()) is None:
                        ch.neighbors[neighbor.get_id()] = {}
                    if len(neighbors[neighbor]):
                        for n in neighbors[neighbor]:
                            if ch.neighbors[neighbor.get_id()].get(n) is None:
                                ch.neighbors[neighbor.get_id()][n.get_id()[1]] = {}
                            ch.neighbors[neighbor.get_id()][n.get_id()[1]] = residue.get_id()[1]
        print('Chain Neighbors Found Successfully')
        return chain.neighbors

    def get_residue_neighbors_shallow(self, residue):
        atoms = self.pdb.get_atoms(residue)
        atoms = self.pdb.probe.get_atoms_in_atom_probe(atoms, residue)
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

    def critical_residues(self, threshold_high, model, chain):
        if not hasattr(self.pdb.structure, 'sasa'):
            self.sasa()
        item = self.pdb.get_item(model, chain)
        if item is None:
            print('----------\nError Getting Chain Neighbors :\nChain not Found\n')
            return None
        critical_residues = self.get_critical_residues(threshold_high, item)
        self.report_critical_residues(critical_residues, chain)

    def get_critical_residues(self, threshold_high, item):
        print('----------\nSearching For Critical Residues', end='\r')
        critical_residues = []
        residues = item.get_residues()
        for residue in residues:
            chain = residue.get_parent()
            model = chain.get_parent()
            if residue.accessibility > float(threshold_high):
                continue
            print('Checking Residue #%s Neighbors' % residue.id[1], end='\r')
            neighbors = self.get_residue_neighbors_shallow(residue)
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
        print('Critical Residues Found Successfully')
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

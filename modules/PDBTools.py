import numpy as np
from modules import FileManager
from models import PDB, Model, Residue, Chain


class PDBTools:

    def __init__(self, pdb_file, probe_points, probe_radius, extended: bool = False, minimal: bool = True):
        self.fm = FileManager(pdb_file, probe_points, probe_radius, extended, minimal)
        self.pdb = self.fm.pdb
        self.fm.ready()

    def sasa(self):
        self.calc_sasa()
        self.fm.add(['sasa'])

    def calc_sasa(self):
        print('----------\nBegin SASA Calculation', end='\r')
        atoms = self.pdb.get_atoms()
        self.pdb.probe.get_neighbor_probe_points(atoms)
        for index, atom in enumerate(atoms):
            pp = self.pdb.probe.points
            atom_probe_points = self.pdb.probe.atoms[index].probe.buried
            atom.accessibility = sum([not p for p in atom_probe_points]) / pp
            atom.size = 4 * np.pi * (atom.radius + self.pdb.probe.radius) ** 2
            atom.sasa = atom.accessibility * atom.size
        self.sum_sasa(self.pdb.structure)
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

    def residue_neighbors(self, model: Model, chain: Chain, residue: Residue):
        item = self.pdb.get_item(model, chain, residue)
        if item is None:
            print('----------\nError Getting Residue Neighbors :\nResidue not Found\n')
            return None
        self.get_residue_neighbors(item)
        self.fm.add(['neighbors'], item)

    def get_residue_neighbors(self, residue):
        print('----------\nSearching Residue Neighbors', end='\r')
        atoms = self.pdb.get_atoms(residue)
        atoms = self.pdb.probe.get_atoms_in_atom_probe(atoms, residue)
        atoms = self.pdb.probe.get_atoms_points_from_neighbors_atoms_probe(atoms, residue)
        residue.neighbors = {'chains': {}}
        for atom in atoms:
            neighbor_residues = {}
            for n in atom.neighbors:
                if neighbor_residues.get(n['atom'].get_parent()) is None:
                    neighbor_residues[n['atom'].get_parent()] = 0
                neighbor_residues.update({n['atom'].get_parent(): neighbor_residues[n['atom'].get_parent()] + n['share']})
            for neighbor_residue in neighbor_residues:
                chain = neighbor_residue.get_parent()
                if residue.neighbors['chains'].get(chain.get_id()) is None:
                    residue.neighbors['chains'][chain.get_id()] = {'residues': {}}
                if residue.neighbors['chains'][chain.get_id()]['residues'].get(neighbor_residue.get_id()[1]) is None:
                    residue.neighbors['chains'][chain.get_id()]['residues'][neighbor_residue.get_id()[1]] = {
                        'name': neighbor_residue.get_resname(),
                        'item': neighbor_residue,
                        'share': 0
                    }
                residue.neighbors['chains'][chain.get_id()]['residues'][neighbor_residue.get_id()[1]]['share'] += neighbor_residues[neighbor_residue]
        print('Residue Neighbors Found Successfully')
        return residue.neighbors

    def chain_neighbors(self, model: Model, chain: Chain):
        item = self.pdb.get_item(model, chain)
        if item is None:
            print('----------\nError Getting Chain Neighbors :\nChain not Found\n')
            return None
        self.get_chain_neighbors(item)
        self.fm.add(['neighbors'], item)

    def get_chain_neighbors(self, chain):
        print('----------\nSearching Chain Neighbors', end='\r')
        chain.neighbors = {'chains': {}}
        for residue in chain.get_residues():
            print('Getting Residue %s%s Neighbors' % (residue.get_resname(), residue.id[1]), end='\r')
            neighbors = self.get_residue_neighbors_shallow(residue)
            if neighbors is not None and len(neighbors['chains']) > 1:
                ch = residue.get_parent()
                for neighbor in [key for key in neighbors['chains'].keys() if key is not ch.get_id()]:
                    if ch.neighbors['chains'].get(neighbor) is None:
                        ch.neighbors['chains'][neighbor] = []
                    if len(neighbors['chains'][neighbor]):
                        for n in neighbors['chains'][neighbor]['residues']:
                            nx = neighbors['chains'][neighbor]['residues'][n]
                            ch.neighbors['chains'][neighbor].append({
                                'from': {residue.get_id()[1]: {'name': residue.get_resname()}},
                                'to': {n: {'name': nx['name']}},
                                'points': nx['points']
                            })
        print('Chain Neighbors Found Successfully')
        return chain.neighbors

    def get_residue_neighbors_shallow(self, residue):
        atoms = self.pdb.get_atoms(residue)
        atoms = self.pdb.probe.get_atoms_in_atom_probe(atoms, residue)
        residue.neighbors = {'chains': {}}
        for atom in atoms:
            neighbor_residues = {}
            for neighbor in atom.neighbors:
                if neighbor_residues.get(neighbor['atom'].get_parent()) is None:
                    neighbor_residues[neighbor['atom'].get_parent()] = 0
                neighbor_residues[neighbor['atom'].get_parent()] += neighbor['share']
            for neighbor_residue in neighbor_residues:
                chain = neighbor_residue.get_parent().get_id()
                if residue.neighbors['chains'].get(chain) is None:
                    residue.neighbors['chains'][chain] = {'residues': {}}
                if residue.neighbors['chains'][chain]['residues'].get(neighbor_residue.get_id()[1]) is None:
                    residue.neighbors['chains'][chain]['residues'][neighbor_residue.get_id()[1]] = {
                        'name': neighbor_residue.get_resname(),
                        'item': neighbor_residue,
                        'points': 0
                    }
                residue.neighbors['chains'][chain]['residues'][neighbor_residue.get_id()[1]]['points'] += neighbor_residues[neighbor_residue]
        return residue.neighbors

    def critical_residues(self, threshold, model:Model, chain:Chain):
        if not hasattr(self.pdb.structure, 'sasa'):
            self.sasa()
        item = self.pdb.get_item(model, chain)
        if item is None:
            print('----------\nError Getting Chain Neighbors :\nChain not Found\n')
            return None
        self.get_critical_residues(threshold, item)
        self.fm.add(['critical'], item)

    def get_critical_residues(self, threshold, chain):
        print('----------\nSearching For Critical Residues', end='\r')
        critical_residues = {}
        for residue in chain.get_residues():
            if residue.accessibility > float(threshold):
                continue
            print('Checking Residue #%s Neighbors' % residue.id[1], end='\r')
            neighbors = self.get_residue_neighbors_shallow(residue)
            siblings = range(residue.get_id()[1] - 3, residue.get_id()[1] + 3)
            equal, non = {}, {}
            for res in neighbors['chains'][chain.get_id()]['residues']:
                if res in siblings:
                    continue
                points = neighbors['chains'][chain.get_id()]['residues'][res]['points']
                if points < 50:
                    continue
                res = neighbors['chains'][chain.get_id()]['residues'][res]['item']
                if res.polar == residue.polar:
                    equal[res.get_id()[1]] = {'name': res.get_resname(), 'points': points}
                else:
                    non[res.get_id()[1]] = {'name': res.get_resname(), 'points': points}
            if len(equal) > 0 or len(non) > 0:
                if residue.polar:
                    if len(equal) == 0:
                        critical_residues[residue.get_id()[1]] = {
                            'name': residue.get_resname(),
                            'sasa': residue.sasa,
                            'type': 'HydPhl-HydPhb',
                            'residues': non
                        }
                    else:
                        critical_residues[residue.get_id()[1]] = {
                            'name': residue.get_resname(),
                            'sasa': residue.sasa,
                            'type': 'HydPhl-HydPhl',
                            'residues': equal
                        }
                else:
                    if len(non) == 0:
                        critical_residues[residue.get_id()[1]] = {
                            'name': residue.get_resname(),
                            'sasa': residue.sasa,
                            'type': 'HydPhb-HydPhb',
                            'residues': equal
                        }

        print('Critical Residues Found Successfully')
        chain.critical = critical_residues
        return chain.critical


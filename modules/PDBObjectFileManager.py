import os
import json
import pickle
from models import PDB, Chain, Model, Residue


class FileManager:
    def __init__(self, pdb_file_address, probe_points, probe_radius, extended: bool = False, minimal: bool = True):
        self.minimal = minimal
        self.extended = extended
        self.pdb = self.__get_pdb(pdb_file_address, probe_points, probe_radius)

    def __get_pdb(self, pdb_file_address, probe_points, probe_radius):
        self.__check_main_folder_existence()
        file, name = self.__get_pdb_file_from_address(pdb_file_address)
        pdb = self.__load_pdb_object(name)
        if pdb is None:
            new = PDB(name, file, probe_points, probe_radius)
            self.__write_pdb_object_to_file(new)
        else:
            new = pdb
        return new

    @staticmethod
    def __get_pdb_file_from_address(address):
        path, file = os.path.split(address)
        name, extension = os.path.splitext(file)
        return file, name

    @staticmethod
    def __check_main_folder_existence():
        if not os.path.isdir(os.getcwd() + '/PDBObject/'):
            os.makedirs(os.getcwd() + '/PDBObject/')

    @staticmethod
    def __load_pdb_object(file):
        if os.path.isfile(os.getcwd() + '/PDBObject/' + file + '.pdbo'):
            with open(os.getcwd() + '/PDBObject/' + file + '.pdbo', 'rb') as f:
                return pickle.load(f)
        else:
            return None

    @staticmethod
    def __duplicate_object_attributes(source, target):
        for k in source.__dict__.keys():
            setattr(target, k, getattr(source, k))

    @staticmethod
    def __write_pdb_object_to_file(pdb):
        with open(os.getcwd() + '/PDBObject/' + pdb.id + '.pdbo', 'wb') as f:
            pickle.dump(pdb, f)

    def ready(self):
        self.__json(self.pdb)
        self.__open()

    def add(self, fields, target=None):
        if target is None:
            target = self.pdb.structure
        for field in fields:
            self.__attach_attribute(field, target)

    def __write(self):
        if self.minimal:
            self.structure = self.__strip_empties_from_dict(self.structure)
        with open(os.getcwd() + './PDBObject/' + self.pdb.structure.id + '.json', 'w') as f:
            f.write(json.dumps(self.structure, default=self.__convert_set_to_list))

    def __open(self) -> bool:
        try:
            file = os.getcwd() + '/PDBObject/' + self.pdb.structure.id + '.json'
            if os.stat(file).st_size == 0:
                return False
            with open(file, 'r') as f:
                self.__append_json(json.load(f))
                return True
        except FileNotFoundError:
            return False

    def __append_json(self, json_structure):
        structure = list(json_structure.keys())[0]
        s = json_structure[structure]
        my_s = self.structure[structure]
        self.__duplicate_fields(my_s, s)
        if 'models' in s:
            for model in s['models'].keys():
                m = s['models'][model]
                my_m = my_s['models'][model]
                self.__duplicate_fields(my_m, m)
                if 'chains' in m:
                    for chain in m['chains'].keys():
                        c = m['chains'][chain]
                        my_c = my_m['chains'][chain]
                        self.__duplicate_fields(my_c, c)
                        if 'residues' in c:
                            for residue in c['residues'].keys():
                                r = c['residues'][residue]
                                my_r = my_c['residues'][residue]
                                self.__duplicate_fields(my_r, r)
                                if 'atoms' in r and self.extended:
                                    for atom in r['atoms'].keys():
                                        a = r['atoms'][atom]
                                        my_a = my_r['atoms'][atom]
                                        self.__duplicate_fields(my_a, a)

    @staticmethod
    def __duplicate_fields(main, target):
        for key in target.keys():
            if key not in main:
                main[key] = target[key]

    def __json(self, pdb):
        self.__structure_to_json(pdb.structure)

    def __add_atoms_to_residue(self, structure):
        if self.extended:
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if 'atoms' not in residue:
                            r = PDB.get_item(Model(model), Chain(chain), Residue(residue))
                            atoms = {}
                            for atom in r:
                                atoms[str(atom.id)] = {}
                            structure[model][chain][residue]['atoms'] = atoms

    def __structure_to_json(self, structure):
        models = {}
        for model in structure:
            chains = {}
            for chain in model:
                residues = {}
                for residue in chain:
                    if self.extended:
                        atoms = {}
                        for atom in residue:
                            atoms[str(atom.id)] = {}
                        residues[str(residue.get_id()[1])] = {'name': residue.get_resname(), 'atoms': atoms}
                    else:
                        residues[str(residue.get_id()[1])] = {'name': residue.get_resname()}
                chains[str(chain.id)] = {'residues': residues}
            models[str(model.id)] = {'chains': chains}
        self.structure = {structure.id: {'models': models}}

    def __attach_attribute(self, field, target):
        level = target.__dict__['level']
        if level == 'S':
            self.__add_field_to_structure(field, target)
        elif level == 'C':
            self.__add_field_to_chain(field, target)
        elif level == 'R':
            self.__add_field_to_residue(field, target)

    def __add_field_to_structure(self, field, target):
        s = self.structure[target.__dict__['_id']]
        for model in target:
            m = s['models'][str(model.id)]
            for chain in model:
                c = m['chains'][str(chain.id)]
                for residue in chain:
                    r = c['residues'][str(residue.get_id()[1])]
                    self.__add_atom_fields_to_residue(residue, r, field)
                    if field in residue.__dict__:
                        r[field] = residue.__dict__[field]
                if field in chain.__dict__:
                    c[field] = chain.__dict__[field]
            if field in model.__dict__:
                m[field] = model.__dict__[field]
        if field in target.__dict__:
            s[field] = target.__dict__[field]
        self.__write()

    def __add_field_to_chain(self, field, target):
        id = target.__dict__['full_id']
        s = self.structure[id[0]]
        m = s['models'][str(id[1])]
        c = m['chains'][str(id[2])]
        if field in target.__dict__:
            c[field] = target.__dict__[field]
        self.__write()

    def __add_field_to_residue(self, field, target):
        id = target.__dict__['full_id']
        s = self.structure[id[0]]
        m = s['models'][str(id[1])]
        c = m['chains'][str(id[2])]
        r = c['residues'][str(id[3][1])]
        self.__add_atom_fields_to_residue(target, r, field)
        if field in target.__dict__:
            r[field] = target.__dict__[field]
        self.__write()

    def __add_atom_fields_to_residue(self, residue, r, field):
        if self.extended:
            for atom in residue:
                a = r['atoms'][str(atom.id)]
                if field in atom.__dict__:
                    a[field] = atom.__dict__[field]

    def __strip_empties_from_list(self, data):
        new_data = []
        for v in data:
            if isinstance(v, dict):
                v = self.__strip_empties_from_dict(v)
            elif isinstance(v, list):
                v = self.__strip_empties_from_list(v)
            if v not in (None, str(), list(), dict(),):
                new_data.append(v)
        return new_data

    def __strip_empties_from_dict(self, data):
        new_data = {}
        for k, v in data.items():
            if isinstance(v, dict):
                if k == 'neighbors' or k == 'critical':
                    new_data[k] = v
                    continue
                v = self.__strip_empties_from_dict(v)
            elif isinstance(v, list):
                if k == 'neighbors' or k == 'critical':
                    new_data[k] = v
                    continue
                v = self.__strip_empties_from_list(v)
            if k == 'name' and len(data.items()) == 1:
                continue
            if v not in (None, str(), list(), dict(),):
                new_data[k] = v
        return new_data

    @staticmethod
    def __convert_set_to_list(obj):
        if isinstance(obj, set):
            return list(obj)
        if 'level' in obj.__dict__:
            level = obj.__dict__['level']
            if level in ['A', 'C', 'M', 'S']:
                return str(obj.get_id())
            elif level == 'R':
                return str(obj.get_id()[1])
        raise TypeError

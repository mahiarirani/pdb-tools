import os
import json


class FileManager:
    def __init__(self, pdb):
        self.pdb = pdb

    def ready(self):
        if not self.__open():
            self.save()

    def save(self):
        self.__json(self.pdb)

    def add(self, fields, target=None):
        if target is None:
            target = self.pdb.structure
        for field in fields:
            self.__attach_attribute(field, target)

    def __write(self):
        with open(os.getcwd() + '/PDBObject/' + self.pdb.structure.id + '.json', 'w') as f:
            f.write(json.dumps(self.structure, default=self.__convert_set_to_list))

    def __open(self) -> bool:
        try:
            file = os.getcwd() + '/PDBObject/' + self.pdb.structure.id + '.json'
            if os.stat(file).st_size == 0:
                return False
            with open(file, 'r') as f:
                self.structure = json.load(f)
                return True
        except FileNotFoundError:
            return False

    def __json(self, pdb):
        self.__structure_to_json(pdb.structure)

    def __structure_to_json(self, structure):
        models = {}
        for model in structure:
            chains = {}
            for chain in model:
                residues = {}
                for residue in chain:
                    atoms = {}
                    for atom in residue:
                        atoms[str(atom.id)] = {}
                    residues[str(residue.get_id()[1])] = {'atoms': atoms}
                chains[str(chain.id)] = {'residues': residues}
            models[str(model.id)] = {'chains': chains}
        self.structure = {structure.id: {'models': models}}
        self.__write()

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
                    for atom in residue:
                        a = r['atoms'][str(atom.id)]
                        if field in atom.__dict__:
                            a[field] = atom.__dict__[field]
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
        for atom in target:
            a = r['atoms'][str(atom.id)]
            if field in atom.__dict__:
                a[field] = atom.__dict__[field]
        if field in target.__dict__:
            r[field] = target.__dict__[field]
        self.__write()

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

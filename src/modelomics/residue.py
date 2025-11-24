# created by clay 11/22/25
'''
the residue class for various protein design tasks
'''

class Residue:
    def __init__(self, atoms, missing_info=None):
        ''' 
        initialize a residue from a list of atoms or missing residue info
        '''
        if atoms:
            # normal residue with atoms
            self.atoms = {atom.name: atom for atom in atoms}
            self._atom_list = list(self.atoms.values())
            first_atom = list(self.atoms.values())[0]
            self.name = first_atom.residue
            self.number = first_atom.resnum
            self.chain = first_atom.chain
            self.insertion = first_atom.insertion
            self.is_missing = False
        elif missing_info:
            # missing residue
            self.atoms = {}
            self._atom_list = []
            self.name = missing_info['name']
            self.number = missing_info['number']
            self.chain = missing_info.get('chain', 'X')
            self.insertion = missing_info.get('insertion', '')
            self.is_missing = True
        else:
            raise ValueError("Either atoms or missing_info must be provided")
        
    def __str__(self):
        if self.is_missing:
            return (f"MISSING RESIDUE {self.name} "
                   f"{self.number}{self.insertion} in chain {self.chain}")
        else:
            return "".join(atom.line for atom in self._atom_list)
    
    def __getitem__(self, key):
        """
        access atom by name (str) or by index (int)
        """
        if self.is_missing:
            raise KeyError(f"Cannot access atoms in missing residue "
                           f"{self.name}{self.number}")
            
        if isinstance(key, int):
            # integer index -> return atom by order
            return self._atom_list[key]
        elif isinstance(key, str):
            # string name -> return atom by dictionary
            return self.atoms[key]
        else:
            raise TypeError(f"Invalid Residue key type: {type(key)}")
        
    def __len__(self):
        return len(self.atoms)

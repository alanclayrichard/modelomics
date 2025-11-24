# created by clay 11/22/25
'''
the residue class for various protein design tasks
'''

class Residue:
    def __init__(self, atoms):
        ''' 
        initialize a residue from a list of atoms 
        '''
        self.atoms = {atom.name: atom for atom in atoms}
        self._atom_list = list(self.atoms.values())
        first_atom = list(self.atoms.values())[0]
        self.name = first_atom.residue
        self.number = first_atom.resnum
        self.chain = first_atom.chain
        
    def __str__(self):
        return (f"Residue {self.name}{self.number} in chain {self.chain} "
                f"with {len(self.atoms)} atoms")
    
    def __getitem__(self, key):
        """
        access atom by name (str) or by index (int)
        """
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

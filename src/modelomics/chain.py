# created by clay 11/22/25
'''
the chain class for use in protein tasks. a chain is a unique molecule
'''

import residue as r

class Chain:
    def __init__(self, atoms):
        ''' 
        initialize the chain object from a list of atoms
        '''
        self.name = atoms[0].chain
        self.atoms = atoms
        self.residues = {}
        self.populate()
        self._residue_list = list(self.residues.values())

    def populate(self):
        """
        efficiently group atoms into residues
        """
            
        resnum = None
        resname = None
        insertion = None
        atoms = []
        
        for atom in self.atoms:
            # check if we're starting a new residue
            if (atom.resnum != resnum or 
                atom.residue != resname or 
                atom.insertion != insertion):
                
                # save previous residue (skip if first)
                if atoms:
                    self.residues[(resname, resnum, insertion)] = r.Residue(
                        atoms
                    )
                
                # start new residue
                resnum = atom.resnum
                resname = atom.residue
                insertion = atom.insertion
                atoms = []
            
            atoms.append(atom)
        
        # add the last residue
        if atoms:
            self.residues[(resname, resnum, insertion)] = r.Residue(atoms)

    def __str__(self):
        return (f"Chain {self.name} with {len(self.residues)} residues "
                f"and {len(self.atoms)} atoms")
    
    def __getitem__(self, key):
        """
        Access residues by:
        - tuple (name, number, insertion) or (name, number)
        - integer index
        """
        if isinstance(key, tuple):
            # handle optional insertion code
            if len(key) == 2:
                key = (key[0], key[1], "")
            return self.residues[key]
        elif isinstance(key, int):
            return self._residue_list[key]
        else:
            raise TypeError(f"Invalid key type: {type(key)}")

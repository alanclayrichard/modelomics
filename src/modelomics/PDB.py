# created by clay 07/01/25
'''
the PDB class for use with protein tools
'''
import atom as a
import chain as c

class PDB:
    def __init__(self, filename):
        '''
        intialize a PDB object from a file
        '''
        self.filename = filename
        self.file_data = ""
        self.atoms = []
        self.chains = {}
        self._chain_list = []
        
        self.parse_lines(filename)
        
    def __str__(self):
        return f"PDB object from{self.filename} with {len(self.chains)} chains"
    
    def __getitem__(self, key):
        '''
        access the chain with a name if str/char and index if int
        '''
        if isinstance(key, str):
            return self.chains[key]
        elif isinstance(key, int):
            return self._chain_list[key]
        else:
            raise TypeError(f"Invalid PDB key type: {type(key)}")
        
    def __len__(self):
        return len(self._chain_list)

    def parse_lines(self, filename):
        '''
        store the atom objects and other data created from lines in the pdb file
        '''
        with open(filename) as f:
            lines = f.readlines()
            
        # Pre-filter lines and create atoms in one pass
        atom_lines = [line for line in lines if line.startswith(("ATOM", "HETATM"))]
        self.atoms = [a.Atom(line) for line in atom_lines]
        
        # Store non-atom lines
        self.file_data = ''.join(line for line in lines if not line.startswith(("ATOM", "HETATM")))
        
        # Build chains efficiently
        self._build_chains()

    def _build_chains(self):
        """
        build chains efficiently in one pass
        """
        self.chains = {}
        current_chain_id = None
        current_atoms = []
        
        for atom in self.atoms:
            if atom.chain != current_chain_id:
                # Save previous chain
                if current_atoms:
                    self.chains[current_chain_id] = c.Chain(current_atoms)
                # Start new chain
                current_chain_id = atom.chain
                current_atoms = []
            current_atoms.append(atom)
        
        # Add last chain
        if current_atoms:
            self.chains[current_chain_id] = c.Chain(current_atoms)

        self._chain_list = list(self.chains.values())

    def strip_waters(self):
        '''
        remove waters from the list of atoms
        '''
        # Filter atoms and rebuild chains in one operation
        self.atoms = [atom for atom in self.atoms if atom.residue != "HOH"]
        self._build_chains()

    def strip_hetatm(self):
        '''
        remove hetatms from the list of atoms
        '''
        # Filter atoms and rebuild chains in one operation  
        self.atoms = [atom for atom in self.atoms if atom.type == "ATOM"]
        self._build_chains()


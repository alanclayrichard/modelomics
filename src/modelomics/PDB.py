# created by clay 07/01/25
'''
the PDB class for use with protein tools
'''
from . import atom as a
from . import chain as c

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
        self.missing_residues = {}
        
        self.parse_lines(filename)
        
    def __str__(self):
        '''
        a descriptive string of the pdb object
        '''
        return f"PDB object from {self.filename} with {len(self.chains)} chains"
    
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
        '''
        the number of chains in this protein
        '''
        return len(self._chain_list)

    def parse_lines(self, filename):
        '''
        store the atom objects and other data created from lines in the pdb file
        '''
        with open(filename) as f:
            lines = f.readlines()
        
        # Parse missing residues - add this block
        for line in lines:
            if line.startswith("REMARK 465") and len(line.split()) >= 5:
                parts = line.split()
                if parts[2] not in ["MISSING", "THE", "FOLLOWING"]:
                    try:
                        res_name = parts[2]
                        chain_id = parts[3]
                        res_num = int(parts[4])
                        if chain_id not in self.missing_residues:
                            self.missing_residues[chain_id] = []
                        self.missing_residues[chain_id].append((res_name, 
                                                                res_num))
                    except (ValueError, IndexError):
                        continue
            
        # pre-filter lines and create atoms in one pass
        atom_lines = [line for line in lines if line.startswith(("ATOM", 
                                                                 "HETATM"))]
        self.atoms = [a.Atom(line) for line in atom_lines]
        
        # store non-atom lines
        self.file_data = ''.join(line for line in lines if not line.startswith((
            "ATOM", "HETATM"
        )))
        
        # build chains efficiently
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
                # save previous chain
                if current_atoms:
                    missing = self.missing_residues.get(current_chain_id, [])
                    self.chains[current_chain_id] = c.Chain(current_atoms, 
                                                            missing)
                # start new chain
                current_chain_id = atom.chain
                current_atoms = []
            current_atoms.append(atom)
        
        # add last chain
        if current_atoms:
            missing = self.missing_residues.get(current_chain_id, [])
            self.chains[current_chain_id] = c.Chain(current_atoms, missing)

        self._chain_list = list(self.chains.values())

    def strip_waters(self):
        '''
        remove waters from the list of atoms
        '''
        # filter atoms and rebuild chains in one operation
        self.atoms = [atom for atom in self.atoms if atom.residue != "HOH"]
        self._build_chains()

    def strip_hetatms(self):
        '''
        remove hetatms from the list of atoms
        '''
        # filter atoms and rebuild chains in one operation  
        self.atoms = [atom for atom in self.atoms if atom.type == "ATOM"]
        self._build_chains()

    def strip_hydrogens(self):
        '''
        remove hydrogen atoms from the list of atoms
        '''
        self.atoms = [atom for atom in self.atoms if atom.element != "H"]
        self._build_chains()

    def rename_to_charmm(self):
        '''
        rename all atoms to charmm conventions
        '''
        # find last residue in each chain for terminal oxygen naming
        last_residues = {}
        for chain_id, chain in self.chains.items():
            if chain._residue_list:
                last_res_num = max(res.number for res in chain._residue_list)
                last_residues[chain_id] = last_res_num
    
        # rename all atoms
        for atom in self.atoms:
            is_last = atom.resnum == last_residues.get(atom.chain, -1)
            atom.rename_to_charmm(is_last)
    
        # rebuild chains to update with new names
        self._build_chains()


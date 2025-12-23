# created by clay 12/19/25
'''
the CIF class for use with protein tools
'''
from . import atom as a
from . import chain as c

class CIF:
    def __init__(self, filename):
        '''
        initialize a CIF object from a mmCIF file
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
        a descriptive string of the CIF object
        '''
        return f"CIF object from {self.filename} with {len(self.chains)} chains"
    
    def __getitem__(self, key):
        '''
        access the chain with a name if str/char and index if int
        '''
        if isinstance(key, str):
            return self.chains[key]
        elif isinstance(key, int):
            return self._chain_list[key]
        else:
            raise TypeError(f"Invalid CIF key type: {type(key)}")
        
    def __len__(self):
        '''
        the number of chains in this protein
        '''
        return len(self._chain_list)

    def parse_lines(self, filename):
        '''
        store the atom objects and other data created from lines in the CIF file
        '''
        with open(filename) as f:
            lines = f.readlines()
        
        # Pre-filter atom lines
        atom_lines = []
        atom_site_idxs = []
        missing_res_idxs = []
        missing_labels = []
        loop_lines = []
        for i, line in enumerate(lines): 
            if line.startswith(("ATOM", "HETATM")):
                atom_lines.append(line)
            elif line.startswith(("_atom_site.")):
                atom_site_idxs.append(line.replace(" ", ""))
                self.file_data += line
            elif line.startswith("_pdbx_unobs_or_zero_occ_residues."):
                missing_res_idxs.append(line.replace(" ", ""))
                missing_labels.append(i)
            elif "loop_" in line:
                loop_lines.append(i)
            else:
                self.file_data += line

        if len(missing_labels) > 0:
            missing_lines = lines[
                missing_labels[-1]+1:
                min([x for x in loop_lines if x > max(missing_labels)])
            ]
            missing_lines = [line for line in missing_lines if '#' not in line]

            missing_res_idxs_dict = {}
            for i, label in enumerate(missing_res_idxs):
                missing_res_idxs_dict[
                    label.replace("_pdbx_unobs_or_zero_occ_residues.", "")[:-1]
                ] = i
            
            for line in missing_lines:
                values = line.split()
                chain_id = values[missing_res_idxs_dict["auth_asym_id"]]
                if chain_id not in self.missing_residues:
                            self.missing_residues[chain_id] = []
                self.missing_residues[chain_id].append(
                    (values[missing_res_idxs_dict["auth_comp_id"]],
                    values[missing_res_idxs_dict["auth_seq_id"]])
                )

        atom_site_idx_dict = {}
        for i, label in enumerate(atom_site_idxs):
            atom_site_idx_dict[label.replace("_atom_site.","")[:-1]] = i
        self.atoms = [a.Atom(line, atom_site_idx_dict) for line in atom_lines]
                
        # build chains efficiently
        self._build_chains()

    def _build_chains(self):
        """
        build chains efficiently in one pass
        """
        self.chains = {}
        chain_atoms = {}
        
        # group all atoms by chain id first
        for atom in self.atoms:
            if atom.chain not in chain_atoms:
                chain_atoms[atom.chain] = []
            chain_atoms[atom.chain].append(atom)
            
        # create chain objects
        for chain_id, atoms in chain_atoms.items():
            missing = self.missing_residues.get(chain_id, [])
            self.chains[chain_id] = c.Chain(atoms, missing)

        self._chain_list = list(self.chains.values())

    def strip_waters(self):
        '''
        remove waters from the list of atoms
        '''
        self.atoms = [atom for atom in self.atoms if atom.residue != "HOH"]
        self._build_chains()

    def strip_hetatms(self):
        '''
        remove hetatms from the list of atoms
        '''
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

    def to_string(self, format='pdb'):
        return ''.join(chain.to_string(format) for chain in self)
    
    def write(self, filename):
        ''' 
        write the contents of this protein to a file
        '''
        with open(filename, 'w') as f:
            if filename.endswith("cif"):
                if self[0][0][0].atom_site_dict is not None:
                    f.write("loop_\n")
                    for label in self[0][0].atom_site_dict:
                        f.write("_atom_site."+label+'\n')
                f.write(self.to_string(format='cif'))
            else: 
                f.write(self.to_string(format='pdb'))


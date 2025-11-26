# created by clay 11/22/25
'''
the chain class for use in protein tasks. a chain is a unique molecule
'''

from . import residue as r
from .utils.aa_names import aa_3to1

class Chain:
    def __init__(self, atoms, missing_residues=None):
        ''' 
        initialize the chain object from a list of atoms and missing residues
        '''
        self.name = atoms[0].chain if atoms else 'X'
        self.atoms = atoms
        self.residues = {}
        self.missing_residues = missing_residues or []
        self.populate()
        self._residue_list = list(self.residues.values())

    def populate(self):
        """
        efficiently group atoms into residues and add missing residues
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
        
        # add the last residue from atoms
        if atoms:
            self.residues[(resname, resnum, insertion)] = r.Residue(atoms)
        
        # now add missing residues as empty residues
        for missing in self.missing_residues:
            res_name = missing[0]
            res_num = missing[1]
            insertion = ""
            
            key = (res_name, res_num, insertion)
            if key not in self.residues:
                self.residues[key] = r.Residue([], missing_info={
                    'name': res_name,
                    'number': res_num,
                    'chain': self.name,
                    'insertion': insertion
                })

    def __str__(self):
        total_residues = len(self.residues)
        missing_count = len(self.missing_residues)
        present_count = total_residues - missing_count
        return (f"Chain {self.name} with {total_residues} residues "
                f"({present_count} present, {missing_count} missing) "
                f"and {len(self.atoms)} atoms")
    
    def __getitem__(self, key):
        """
        access residues by:
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
        
    @property
    def sequence(self):
        """
        return full sequence including missing residues, ordered by number
        """
        # sort residues by number to get proper sequence order
        sorted_residues = sorted(self.residues.items(), 
                                 key=lambda x: x[0][1])
        return "".join(aa_3to1.get(res_info[0], "X") 
                       for res_info, res_obj in sorted_residues)
    
    def write_pdb(self, filename):
        ''' 
        write the contents of this protein to a pdb file
        '''
        with open(filename, 'w') as f:
            for atom in self.atoms:
                f.write(atom.line)

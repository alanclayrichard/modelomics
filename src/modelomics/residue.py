# created by clay 11/22/25
'''
the residue class for various protein design tasks
'''

import math

from . import atom as a
from .vdw_point_cloud import VDWPointCloud

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
        
    @property
    def pos(self):
        """
        get xyz coordinates of alpha carbon or average of all atoms
        """
        if self.is_missing:
            return [0.0, 0.0, 0.0]
            
        if 'CA' in self.atoms:
            # return alpha carbon position
            ca = self.atoms['CA']
            return [ca.x, ca.y, ca.z]
        elif self.atoms:
            # return average position of all atoms
            x_sum = sum(atom.x for atom in self._atom_list)
            y_sum = sum(atom.y for atom in self._atom_list)
            z_sum = sum(atom.z for atom in self._atom_list)
            n_atoms = len(self._atom_list)
            return [x_sum/n_atoms, y_sum/n_atoms, z_sum/n_atoms]
        else:
            return [0.0, 0.0, 0.0]
    
    def covalent_bond(self, other_residue, tolerance=0.2):
        """
        determine whether this residue is covalently bonded to another residue
        using atom-specific covalent radii.
        """
        if self.is_missing or other_residue.is_missing:
            return False

        for atom1 in self._atom_list:
            r1 = atom1.radius
            for atom2 in other_residue._atom_list:
                r2 = atom2.radius

                # expected maximum covalent-bond distance
                max_dist = r1 + r2

                dx = atom1.x - atom2.x
                dy = atom1.y - atom2.y
                dz = atom1.z - atom2.z
                dist = (dx*dx + dy*dy + dz*dz) ** 0.5

                if dist <= 0.6*max_dist:
                    return True

        return False

    
    def hydrogen_bond(self, other_residue, ha_cutoff=2.5, min_dha_angle=120):
        """
        determine if this residue has a hydrogen bond to another residue
        criteria: H...A distance < 3.0A, D-H...A angle > 120deg
        """
        if self.is_missing or other_residue.is_missing:
            return False
        
        # simplified: any nitrogen or oxygen can be donor/acceptor
        def is_polar_atom(atom):
            return (hasattr(atom, 'element') and atom.element in ['N', 'O']) or \
                   atom.name.startswith(('N', 'O'))
        
        def angle_between_vectors(v1, v2):
            """calculate angle between two vectors in degrees"""
            dot_product = sum(a*b for a, b in zip(v1, v2))
            mag1 = sum(a**2 for a in v1)**0.5
            mag2 = sum(a**2 for a in v2)**0.5
            if mag1 == 0 or mag2 == 0:
                return 0
            cos_angle = dot_product / (mag1 * mag2)
            cos_angle = max(-1, min(1, cos_angle))
            return math.degrees(math.acos(cos_angle))
        
        # check both directions: self - other and other - self
        for res1, res2 in [(self, other_residue), (other_residue, self)]:
            for donor in res1._atom_list:
                if not is_polar_atom(donor):
                    continue
                    
                # find hydrogen bonded to donor
                h_atom = None
                for atom in res1._atom_list:
                    if atom.element == 'H':
                        h_dist = ((donor.x - atom.x)**2 + 
                                 (donor.y - atom.y)**2 + 
                                 (donor.z - atom.z)**2)**0.5
                        if h_dist < 1.2:
                            h_atom = atom
                            break
                
                if not h_atom:
                    continue
                    
                for acceptor in res2._atom_list:
                    if not is_polar_atom(acceptor):
                        continue
                        
                    # check H...A distance
                    ha_dist = ((h_atom.x - acceptor.x)**2 + 
                              (h_atom.y - acceptor.y)**2 + 
                              (h_atom.z - acceptor.z)**2)**0.5
                                              
                    if ha_dist > ha_cutoff:
                        continue
                        
                    # check D-H...A angle (angle at hydrogen atom)
                    hd_vec = [donor.x - h_atom.x, 
                             donor.y - h_atom.y, 
                             donor.z - h_atom.z]
                    ha_vec = [acceptor.x - h_atom.x, 
                             acceptor.y - h_atom.y, 
                             acceptor.z - h_atom.z]
                             
                    dha_angle = angle_between_vectors(hd_vec, ha_vec)
                    
                    if dha_angle > min_dha_angle:
                        # print(f"hbond found: {res1.name}{res1.number} "
                        #       f"{donor.name}-{h_atom.name}..."
                        #       f"{res2.name}{res2.number} {acceptor.name}, "
                        #       f"ha_dist={ha_dist:.2f}, angle={dha_angle:.1f}")
                        return True
        return False
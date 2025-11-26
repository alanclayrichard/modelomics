# created by clay 07/01/25
'''
the atom class to be used in various protein tasks
'''

from .utils.topology import (vdwRadii, 
                             _POLAR_ELEMENTS, 
                             _POSITIVE_ATOMS, 
                             _NEGATIVE_ATOMS, 
                             _CTERM_NEGATIVE, 
                             _NTERM_POSITIVE
)

class Atom:
    def __init__(self, line):
        self.type = line[:6].strip()
        self.number = int(line[6:11])
        self.name = line[12:16].strip()
        self.residue = line[17:20].strip()
        self.chain = line[21]
        self.resnum = int(line[22:26])
        self.insertion = line[26].replace(" ", "")
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        self.occupancy = float(line[54:60])
        self.temp_factor = float(line[60:66])
        self.element = line[76:78].strip()
        self.radius = vdwRadii[self.element.upper()]
        self.line = line
        
        # compute properties
        self.polar = self._is_polar()
        self.charge = self._get_charge()

    def __str__(self):
        return self.line
    
    def _is_polar(self):
        return self.element.upper() in _POLAR_ELEMENTS
    
    def _get_charge(self):
        """
        determine formal charge based on residue type and atom name
        returns: +1 for positive, -1 for negative, 0 for neutral
        """
        # check positive charges
        if self.residue in _POSITIVE_ATOMS:
            if self.name in _POSITIVE_ATOMS[self.residue]:
                return 1
                
        # check negative charges  
        if self.residue in _NEGATIVE_ATOMS:
            if self.name in _NEGATIVE_ATOMS[self.residue]:
                return -1
                
        # check n-terminus (first residue in chain)
        if self.resnum == 1 and self.name in _NTERM_POSITIVE:
            return 1
            
        # check c-terminus 
        if self.name == "OXT":
            return -1
            
        return 0  # neutral
    
    def rename_to_charmm(self, is_last_residue=False):
        """
        rename atom names to follow CHARMM conventions
        based on common PDB -> CHARMM naming conversions
        """
        original_name = self.name
        
        # isoleucine delta carbon
        if self.residue == "ILE" and self.name == "CD1":
            self.name = "CD"
        
        # terminal oxygens
        elif is_last_residue and self.name == "O":
            self.name = "OT1"
        elif is_last_residue and self.name == "OXT":
            self.name = "OT2"
        
        # n-terminal hydrogens
        elif self.name in ["1H", "2H", "3H"]:
            digit = self.name[0]
            self.name = f"HT{digit}"
        
        # backbone hydrogen
        elif self.name == "H":
            self.name = "HN"
        
        # general hydrogen naming: move leading digit to end
        elif self._is_hydrogen() and self.name[0].isdigit():
            digit = self.name[0]
            rest = self.name[1:]
            self.name = f"{rest}{digit}"
        
        # asparagine and glutamine sidechain oxygens
        elif self.residue == "ASN" and self.name == "OD1":
            self.name = "OD1"  # already correct
        elif self.residue == "GLN" and self.name == "OE1":
            self.name = "OE1"  # already correct
        
        # proline hydrogens
        elif (self.residue == "PRO" and 
              self.name.startswith("H") and 
              self.name[1:].startswith(("G", "D", "B"))):
            # move digit from start to end for proline hydrogens
            if len(self.name) > 2 and self.name[0].isdigit():
                digit = self.name[0]
                rest = self.name[1:]
                self.name = f"{rest}{digit}"
        
        return original_name != self.name  # return True if name was changed

    def _is_hydrogen(self):
        """check if atom is a hydrogen"""
        return self.element.upper() == "H" or self.name.startswith("H")

    def rename_from_charmm(self, is_last_residue=False):
        """
        rename atom names from CHARMM back to PDB conventions
        reverses the rename_to_charmm function
        """
        original_name = self.name
        
        # isoleucine delta carbon
        if self.residue == "ILE" and self.name == "CD":
            self.name = "CD1"
        
        # terminal oxygens
        elif is_last_residue and self.name == "OT1":
            self.name = "O"
        elif is_last_residue and self.name == "OT2":
            self.name = "OXT"
        
        # n-terminal hydrogens
        elif self.name.startswith("HT") and len(self.name) == 3:
            digit = self.name[2]
            if digit.isdigit():
                self.name = f"{digit}H"
        
        # backbone hydrogen
        elif self.name == "HN":
            self.name = "H"
        
        # general hydrogen naming: move trailing digit to start
        elif (self._is_hydrogen() and 
              len(self.name) > 1 and 
              self.name[-1].isdigit()):
            digit = self.name[-1]
            rest = self.name[:-1]
            self.name = f"{digit}{rest}"
        
        return original_name != self.name  # return True if name was changed
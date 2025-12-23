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
    def __init__(self, line, atom_site_dict = None):
        self.parsed_line = line
        self.atom_site_dict = atom_site_dict
        self.from_cif = atom_site_dict is not None
        
        if atom_site_dict != None:
            self._parse_cif_line(line, atom_site_dict)
        else:
            self._parse_pdb_line(line)
        
        # compute properties
        self.polar = self._is_polar()
        self.charge = self._get_charge()
        
        # compute properties
        self.polar = self._is_polar()
        self.charge = self._get_charge()

    def __str__(self):
        return self.to_string()
    
    def _parse_pdb_line(self, line):
        """parse standard PDB ATOM/HETATM line"""
        self.type = line[:6].strip()
        self.number = int(line[6:11])
        self.name = line[12:16].strip()
        self.altloc = line[16]
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
        self.radius = vdwRadii.get(self.element.upper(), 0)

    def _parse_cif_line(self, line, atom_site_dict):
        '''parse cif atom lines given column mapping'''
        tokens = line.split()
        try:
            self.type = tokens[atom_site_dict["group_PDB"]]
            self.number = int(tokens[atom_site_dict["id"]])
            self.name = tokens[atom_site_dict["label_atom_id"]]
            alt_idx = atom_site_dict.get("label_alt_id")
            if alt_idx is not None:
                val = tokens[alt_idx]
                self.altloc = " " if val in (".", "?") else val
            else:
                self.altloc = " "
            self.element = tokens[atom_site_dict["type_symbol"]]
            self.residue = tokens[atom_site_dict["auth_comp_id"]]
            self.chain = tokens[atom_site_dict["auth_asym_id"]]
            label_idx = atom_site_dict.get("auth_seq_id")
            resnum = -1
            if label_idx is not None:
                val = tokens[label_idx]
                if val not in (".", "?"):
                    resnum = int(val)
            self.resnum = resnum
            ins_idx = atom_site_dict.get("pdbx_PDB_ins_code")
            if ins_idx is not None:
                ins = tokens[ins_idx]
                self.insertion = "" if ins in ("?", ".") else ins
            else:
                self.insertion = ""

            self.x = float(tokens[atom_site_dict["Cartn_x"]])
            self.y = float(tokens[atom_site_dict["Cartn_y"]])
            self.z = float(tokens[atom_site_dict["Cartn_z"]])

            occ_idx = atom_site_dict.get("occupancy")
            self.occupancy = float(tokens[occ_idx]) if occ_idx is not None else 1.0

            b_idx = atom_site_dict.get("B_iso_or_equiv")
            self.temp_factor = float(tokens[b_idx]) if b_idx is not None else 0.0

        except (IndexError, ValueError) as e:
            raise ValueError(f"Cannot parse CIF line: {line}") from e

        self.radius = vdwRadii.get(self.element.upper(), 0)

    def to_string(self, format='pdb'):
        if format == 'cif':
            if self.from_cif:
                tokens = self.parsed_line.split()
                d = self.atom_site_dict

                tokens[d["group_PDB"]] = self.type
                tokens[d["id"]] = str(self.number)
                tokens[d["label_atom_id"]] = self.name

                alt_idx = d.get("label_alt_id")
                if alt_idx is not None:
                    tokens[alt_idx] = self.altloc if self.altloc != " " else "."

                tokens[d["type_symbol"]] = self.element
                tokens[d["auth_comp_id"]] = self.residue
                tokens[d["auth_asym_id"]] = self.chain
                tokens[d["auth_seq_id"]] = str(self.resnum) if self.resnum != -1 else "?"
                tokens[d["pdbx_PDB_ins_code"]] = self.insertion if self.insertion else "?"

                tokens[d["Cartn_x"]] = f"{self.x:.3f}"
                tokens[d["Cartn_y"]] = f"{self.y:.3f}"
                tokens[d["Cartn_z"]] = f"{self.z:.3f}"
                tokens[d["occupancy"]] = f"{self.occupancy:.2f}"
                tokens[d["B_iso_or_equiv"]] = f"{self.temp_factor:.2f}"

                return " ".join(tokens) + "\n"

            else:
                self.atom_site_dict = {k: i for i, k in enumerate([
                    "group_PDB",
                    "id",
                    "type_symbol",
                    "label_atom_id",
                    "label_alt_id",
                    "label_comp_id",
                    "label_asym_id",
                    "label_seq_id",
                    "auth_comp_id",
                    "auth_asym_id",
                    "auth_seq_id",
                    "pdbx_PDB_ins_code",
                    "Cartn_x",
                    "Cartn_y",
                    "Cartn_z",
                    "occupancy",
                    "B_iso_or_equiv",
                ])}
                return (
                    f"{self.type} "
                    f"{self.number} "
                    f"{self.element} "
                    f"{self.name} "
                    f"{self.altloc if self.altloc != ' ' else '.'} "
                    f"{self.residue} "        # label_comp_id
                    f"{self.chain} "          # label_asym_id
                    f"{self.resnum} "         # label_seq_id
                    f"{self.residue} "        # auth_comp_id
                    f"{self.chain} "          # auth_asym_id
                    f"{self.resnum} "         # auth_seq_id
                    f"{self.insertion if self.insertion else '?'} "
                    f"{self.x:.3f} "
                    f"{self.y:.3f} "
                    f"{self.z:.3f} "
                    f"{self.occupancy:.2f} "
                    f"{self.temp_factor:.2f}\n"
                )


        elif format == 'pdb':
            return (
                f"{self.type:<6}"          # 1-6
                f"{self.number:>5} "       # 7-11 + 12
                f"{self.name:<4}"          # 13-16
                f"{self.altloc:1}"         # 17 
                f"{self.residue:<4}"       # 18-21 
                f"{self.chain:1}"          # 22
                f"{self.resnum:>4}"        # 23-26
                f"{self.insertion:<1}   "  # 27 + 28-30
                f"{self.x:>8.3f}"          # 31-38
                f"{self.y:>8.3f}"          # 39-46
                f"{self.z:>8.3f}"          # 47-54
                f"{self.occupancy:>6.2f}"  # 55-60
                f"{self.temp_factor:>6.2f}"# 61-66
                f"          "              # 67-76
                f"{self.element:>2}\n"     # 77-78
            )
        else:
            raise ValueError(f"format type {filetype} not known")

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
    
    def move(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

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
    
    def covalent_bond(self, other_atom, tolerance=0.4):
        """
        check for bonding using a standard buffer (tolerance) in Angstroms.
        standard covalent bond distance is roughly sum of radii.
        """
        max_dist = self.radius + other_atom.radius + tolerance
        min_dist = 0.4 
        
        dx = self.x - other_atom.x
        dy = self.y - other_atom.y
        dz = self.z - other_atom.z
        dist_sq = dx**2 + dy**2 + dz**2

        return (min_dist**2) <= dist_sq <= (max_dist**2)
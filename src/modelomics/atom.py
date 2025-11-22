# created by clay 07/01/25
'''
the atom class to be used in various protein tasks
'''

from utils.vdw_radii import vdwRadii

# polar elements for faster lookups
_POLAR_ELEMENTS = {"N", "O", "F", "S", "CL"}

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
        self.radius = vdwRadii[self.element]
        self.line = line

    def __str__(self):
        return self.line
    
    def _is_polar(self):
        return self.element.upper() in _POLAR_ELEMENTS
class Atom:
    def __init__(self, line):
        self.number = int(line[6:11].strip())
        self.name = line[12:16].strip()
        self.residue = line[17:20].strip()
        self.chain = line[21]
        self.resnum = int(line[22:26].strip())
        self.x = float(line[30:38].strip())
        self.y = float(line[38:46].strip())
        self.z = float(line[46:54].strip())
        self.occupancy = float(line[54:60].strip())
        self.temp_factor = float(line[60:66].strip())
        self.element = line[76:78].strip()
        self.line = line

    def __str__(self):
        return self.line

    
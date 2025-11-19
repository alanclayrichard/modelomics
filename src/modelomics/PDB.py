import atom
class PDB:
    def __init__(self, filename):
        self.filename = filename
        self.atoms = self.parse_lines(filename)

    def parse_lines(self, filename):
        atoms = []
        with open(filename) as f:
            lines = f.readlines()
            atoms.extend([atom.Atom(line) for line in lines if line.startswith(("ATOM", "HETATM"))])
        return atoms

    

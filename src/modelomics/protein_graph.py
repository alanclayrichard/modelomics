# created by clay 07/01/25
import os
from Bio.PDB import PDBParser, MMCIFParser

# base graph builder class
class BaseProteinGraphBuilder:
    def __init__(self, cutoff=8.0):
        # store cutoff distance
        self.cutoff = cutoff

    def _parse_structure(self, file_path):
        # get file extension
        ext = os.path.splitext(file_path)[-1].lower()

        # select parser
        if ext == ".pdb":
            parser = PDBParser(QUIET=True)
        elif ext == ".cif":
            parser = MMCIFParser(QUIET=True)
        else:
            raise ValueError(f"Unsupported file type: {ext}")

        # parse structure
        structure = parser.get_structure("structure", file_path)
        return structure

    def build(self, filename, chain):
        # this method should be implemented in subclasses
        raise NotImplementedError("Must be implemented by subclass.")
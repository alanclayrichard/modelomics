
# created by clay 07/01/25
import numpy as np
from scipy.spatial import KDTree

import torch
from torch_geometric.data import Data

from .protein_graph import BaseProteinGraphBuilder

# define atomic radii for common elements
AtomicRadii = {
    "H": 1.0, "C": 1.7, "N": 1.625, "O": 1.48, "S": 1.782
}

# atom-level graph builder
class AtomGraphBuilder(BaseProteinGraphBuilder):
    def __init__(self, cutoff=8.0):
        # call parent constructor
        super().__init__(cutoff)
        # initialize element and residue encoders
        self.element_to_id = {}
        self.residue_to_id = {}

    def _encode_categories(self, items, lookup):
        # encode items to integers using a lookup dictionary
        ids = []
        for item in items:
            if item not in lookup:
                lookup[item] = len(lookup)
            ids.append(lookup[item])
        return np.array(ids), lookup

    def _parse_atoms(self, structure, chain_id):
        # collect atom-level features
        positions, serials, radii = [], [], []
        elements, residues = [], []

        for model in structure:
            for chain in model:
                if chain_id and chain.id != chain_id:
                    continue
                for residue in chain:
                    resname = residue.get_resname()
                    for atom in residue:
                        coord = atom.coord
                        serial = atom.serial_number
                        element = atom.element.strip().capitalize()

                        positions.append(coord)
                        serials.append(serial)
                        radii.append(AtomicRadii.get(element, 0.0))
                        elements.append(element)
                        residues.append(resname)

        return (
            np.array(positions),
            np.array(serials),
            np.array(radii),
            elements,
            residues
        )

    def _build_edges(self, positions):
        # build KD-tree and find neighbor pairs within cutoff
        tree = KDTree(positions)
        pairs = tree.query_pairs(r=self.cutoff, output_type='ndarray')
        edges = np.vstack([pairs, pairs[:, [1, 0]]]).T
        return edges

    def build(self, filename, chain=None):
        # parse structure from file
        structure = self._parse_structure(filename)

        # extract atom-level data
        pos_np, serials_np, radii_np, elements, residues = self._parse_atoms(
            structure, 
            chain
        )

        # encode element and residue names to integers
        element_ids, self.element_to_id = self._encode_categories(
            elements, 
            self.element_to_id
        )
        residue_ids, self.residue_to_id = self._encode_categories(
            residues, 
            self.residue_to_id
        )

        # stack features: serial, radius, element ID, residue ID
        features_np = np.stack([
            serials_np,
            radii_np,
            element_ids,
            residue_ids
        ], axis=1)

        # convert to torch tensors
        x = torch.tensor(features_np, dtype=torch.float)
        pos = torch.tensor(pos_np, dtype=torch.float)
        edge_index = torch.tensor(self._build_edges(pos_np), dtype=torch.long)

        # return a PyG Data object
        return Data(x=x, edge_index=edge_index, pos=pos)
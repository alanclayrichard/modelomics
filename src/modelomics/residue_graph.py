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

# residue-level graph builder
class ResidueGraphBuilder(BaseProteinGraphBuilder):
    def __init__(self, cutoff=8.0):
        # call parent constructor
        super().__init__(cutoff)

    def _parse_residues(self, structure, chain_id):
        # collect residue-level data
        residues = []

        for model in structure:
            for chain in model:
                if chain_id and chain.id != chain_id:
                    continue

                for residue in chain:
                    atom_coords = []
                    atom_elements = []
                    atom_radii = []
                    ca_coord = None

                    for atom in residue:
                        coord = atom.coord
                        element = atom.element.strip().capitalize()
                        radius = AtomicRadii.get(element, 0.0)

                        if atom.get_id() == "CA":
                            ca_coord = coord

                        atom_coords.append(coord)
                        atom_elements.append(element)
                        atom_radii.append(radius)

                    if len(atom_coords) == 0:
                        continue

                    residues.append({
                        "resname": residue.get_resname(),
                        "ca_coords": ca_coord,
                        "atom_coords": np.array(atom_coords),
                        "atom_elements": atom_elements,
                        "atom_radii": atom_radii,
                        "chain": chain.id
                    })

        return residues

    def _build_edges(self, residues):
        # use mean atom coords as residue centroids
        centroids = np.array([r["atom_coords"].mean(axis=0) for r in residues])
        tree = KDTree(centroids)
        pairs = tree.query_pairs(r=self.cutoff, output_type='ndarray')
        edge_index = np.vstack([pairs, pairs[:, [1, 0]]]).T
        return edge_index

    def build(self, filename, chain=None):
        # parse structure from file
        structure = self._parse_structure(filename)

        # extract residue-level info
        residues = self._parse_residues(structure, chain)

        # build graph edges between residues
        edge_index = torch.tensor(self._build_edges(residues), dtype=torch.long)

        # create node features: CA coordinates for each residue
        ca_coords = []
        for r in residues:
            if r["ca_coords"] is not None:
                ca_coords.append(r["ca_coords"])
            else:
                ca_coords.append(np.zeros(3))  # fallback

        x = torch.tensor(np.array(ca_coords), dtype=torch.float)

        # keep detailed residue info if desired
        node_features = []
        for r in residues:
            node_features.append({
                "resname": r["resname"],
                "atom_coords": torch.tensor(r["atom_coords"], dtype=torch.float),
                "atom_elements": r["atom_elements"],
                "atom_radii": torch.tensor(r["atom_radii"], dtype=torch.float),
                "chain": r["chain"]
            })

        # build PyG Data object with x and pos set
        data = Data(x=x, edge_index=edge_index)
        data.pos = x  # for test compatibility
        data.residues = node_features
        return data
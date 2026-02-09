# created by clay 02/09/26
'''
create an all-atom graph representation of a protein chain(s)
'''
# python imports
import torch
from torch_geometric.data import Data
from torch_geometric.nn import radius_graph
import numpy as np

# local imports
from .utils.atom_names import charmm_atom_names, charmm_atom_name_to_idx, CHARMM_ATOM_ONEHOT

class AtomGraph(Data):
    def __init__(self, chains=None, cutoff=5.0):
        if chains is None:
            super().__init__()
            return

        # accept a single chain or a list of chains
        if not isinstance(chains, list):
            chains = [chains]

        x, pos, edge_index, edge_attr = self._construct_graph(chains, cutoff)
        super().__init__(
            x=torch.tensor(x, dtype=torch.float),
            edge_index=torch.tensor(edge_index, dtype=torch.long),
            edge_attr=torch.tensor(edge_attr, dtype=torch.float),
            pos=torch.tensor(pos, dtype=torch.float),
        )
        self.chains = chains

    def _construct_graph(self, chains, cutoff):
        x, pos = self._nodes_pos(chains)
        edge_index, edge_attr = self._edges(pos, cutoff)
        return x, pos, edge_index, edge_attr

    def _nodes_pos(self, chains):
        '''
        construct node features and positions from all atoms in chains.
        node features are one-hot encoded CHARMM atom names.
        '''
        # collect all atoms across chains
        atoms = []
        for chain in chains:
            atoms.extend(chain.atoms)

        node_features = []
        coords = []
        for atom in atoms:
            idx = charmm_atom_name_to_idx.get(atom.name)
            if idx is not None:
                node_features.append(CHARMM_ATOM_ONEHOT[idx])
            else:
                node_features.append(np.zeros(len(charmm_atom_names)))
            coords.append([atom.x, atom.y, atom.z])

        return np.array(node_features), np.array(coords)

    def _edges(self, pos, cutoff):
        '''
        build directed edges within radius cutoff with distance and
        unit vector direction as features
        '''
        pos_t = torch.tensor(pos, dtype=torch.float)

        edge_index = radius_graph(
            pos_t,
            r=cutoff,
            loop=False
        )

        src, dst = edge_index
        diff = pos_t[dst] - pos_t[src]
        dist = torch.norm(diff, dim=1, keepdim=True)
        unit_vec = diff / dist.clamp(min=1e-8)

        # edge_attr: [distance, ux, uy, uz]
        edge_attr = torch.cat([dist, unit_vec], dim=1)

        # sort edges by distance
        order = dist.squeeze().argsort()
        edge_index = edge_index[:, order]
        edge_attr = edge_attr[order]

        return edge_index.numpy(), edge_attr.numpy()

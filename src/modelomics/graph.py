# created by clay 07/01/25
'''
create a graph representation of a protein via the modelomics.pdb object
'''
# python imports
import os
import torch
from torch_geometric.data import Data
from torch_geometric.nn import radius_graph
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=SyntaxWarning)
from prody import parsePDB, GNM, calcSqFlucts

# local imports
from .vdw_point_cloud import VDWPointCloud
from .utils.aa_names import *
from .utils.topology import hydrophobicity

class Graph(Data):
    def __init__(self, chain=None, cutoff=15.0):
        if chain is None:
            super().__init__()
            return

        x, pos, edge_index, edge_attr = self._construct_graph(chain, cutoff)
        super().__init__(
            x=torch.tensor(x, dtype=torch.float),
            edge_index=torch.tensor(edge_index, dtype=torch.long),
            edge_attr=torch.tensor(edge_attr, dtype=torch.float),
            pos=torch.tensor(pos, dtype=torch.float),
        )
        self.chain = chain

    def _construct_graph(self, chain, cutoff):
        x, pos = self._nodes_pos(chain)
        edge_index, edge_attr = self._edges(chain, cutoff)
        return x, pos, edge_index, edge_attr

    def _nodes_pos(self, chain):
        '''
        construct the node features and get the pos tensor
        '''
        node_features = []

        # sequence length
        L = len(chain._residue_list)
        
        # sasa precalculation
        pc = VDWPointCloud(chain.atoms, N=92, probe_radius=1.4)
        pc.populate()
        pc.remove_buried()
        
        # normal mode fluctuations 
        # write tmp file
        chain.write_pdb("./tmp.pdb")
        nma_vals = self._get_nma_values("./tmp.pdb", chain.name)
        # remove tmp pdb
        if os.path.exists("./tmp.pdb"):
            os.remove("./tmp.pdb")

        coords = []            
        for i, res in enumerate(chain._residue_list):
            res_feats = []
            res_1 = convert_3to1(res.name)

            # get the position of this node (CA xyz, or average if no CA)
            ca = res.pos
            coords.append(ca)

            # the one hot encoded aa
            aa_idx = aa_list.index(res_1)
            res_feats.extend(AA_ONEHOT[aa_idx])
            
            # physicochemical properties (binary) 5 dimensions
            res_feats.extend([
                1 if res_1 in polar_aa else 0,
                1 if res_1 in nonpolar_aa else 0,
                1 if res_1 in positive_aa else 0,
                1 if res_1 in negative_aa else 0,
                1 if res_1 in aromatic_aa else 0,
                1 if res_1 in aliphatic_aa else 0
            ])

            # the internal, relative index normalized by length
            res_feats.append(i/L)

            # the hydrophobicity of the aa 
            res_feats.append(hydrophobicity[res_1])
            
            # sasa features 5 dimensions
            sasa = pc.residue_sasa_by_type([res])
            res_feats.extend([
                sasa['total'],
                sasa['polar'],
                sasa['nonpolar'],
                sasa['positive'],
                sasa['negative']
            ])
            
            # nma features 1 dimension
            res_feats.append(nma_vals.get(res.number, 0.0))
            
            node_features.append(res_feats)

        return np.array(node_features), np.array(coords)
    
    def _edges(self, chain, cutoff):
        pos = torch.tensor(
            np.array([res.pos for res in chain._residue_list]),
            dtype=torch.float
        )

        edge_index = radius_graph(
            pos,
            r=cutoff,
            loop=False
        )

        # Edge attributes (scaled distances)
        src, dst = edge_index
        coord_diff = pos[src] - pos[dst]
        dist = torch.norm(coord_diff, dim=1) / cutoff

        # Precompute bond maps once
        covalent = torch.zeros(len(dist))
        hbond = torch.zeros(len(dist))

        for k, (i, j) in enumerate(zip(src.tolist(), dst.tolist())):
            resi = chain._residue_list[i]
            resj = chain._residue_list[j]
            covalent[k] = float(resi.covalent_bond(resj))
            hbond[k] = float(resi.hydrogen_bond(resj, hydrogen_present=False))

        edge_attr = torch.stack([dist, covalent, hbond], dim=1)
        return edge_index.numpy(), edge_attr.numpy()

    
    def _get_nma_values(self, pdb_path, chain_id):
        '''
        get the normal mode MSFs from gaussian network model
        (motion in all directions)
        '''
        try:
            structure = parsePDB(pdb_path)
            if structure is None: return {}
            calphas = structure.select(f'chain {chain_id} and calpha')
            if calphas is None: return {}
            
            gnm = GNM(f'gnm')
            gnm.buildKirchhoff(calphas)
            gnm.calcModes()
            sq_flucts = calcSqFlucts(gnm)
            
            return dict(zip(calphas.getResnums(), sq_flucts))
        except Exception as e:
            print(f"NMA failed for {pdb_path}_{chain_id}: {e}")
            return {}
 
# created by clay 07/01/25
'''
create a graph representation of a protein via the modelomics.pdb object
'''
# python imports
import os
import torch
from torch_geometric.data import Data
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=SyntaxWarning)
from prody import parsePDB, GNM, calcSqFlucts, confProDy

# local imports
from .vdw_point_cloud import VDWPointCloud
from .sequence_embeddings import embed_sequence
from .utils.aa_names import *

class Graph(Data):
    def __init__(self, chain, cutoff=15.0):
        # build nodes, edges, positions, and features
        x, pos, edge_index, edge_attr = self._construct_graph(chain, cutoff)

        # initialize the parent class (torch geometric graph object)
        super().__init__(
            x=torch.tensor(x, dtype=torch.float), 
            edge_index=torch.tensor(edge_index, dtype=torch.long), 
            edge_attr=torch.tensor(edge_attr, dtype=torch.float), 
            pos=torch.tensor(pos, dtype=torch.float)
        )
        
        # protein objects 
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
        
        # sequence embeddings  
        try:
            emb = embed_sequence(chain.sequence)
            if isinstance(emb, torch.Tensor):
                emb = emb.to(torch.float32).detach().cpu().numpy()
            # handle both ESMC (3D: batch, seq+boundary, dim) and E1 (2D: seq, dim)
            if emb.ndim == 3:
                emb = emb.squeeze(0)  # remove batch dimension
                # remove boundary tokens for ESMC (first and last)
                if emb.shape[0] == len(chain.sequence) + 2:
                    emb = emb[1:-1]
        except Exception as e:
            print(f"Embedding failed for {chain} with {e}")
            print(f"N residues: {len(chain)} "
                  f"N aa_seq: {len(chain.sequence)}")
            emb = np.zeros((len(chain.sequence), 960))

        coords = []            
        for i, res in enumerate(chain._residue_list):
            res_feats = []
            res_1 = convert_3to1(res.name)

            # get the position of this node (CA xyz, or average if no CA)
            ca = res.pos
            coords.append(ca)
            
            # embeddings 960 dimensions
            res_feats.extend(emb[i])
            
            # physicochemical properties (binary) 5 dimensions
            res_feats.extend([
                1 if res_1 in polar_aa else 0,
                1 if res_1 in nonpolar_aa else 0,
                1 if res_1 in positive_aa else 0,
                1 if res_1 in negative_aa else 0,
                1 if res_1 in aromatic_aa else 0
            ])
            
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
        # calculate pairwise distances
        edge_index = []
        edge_attr = []
        
        for i, resi in enumerate(chain._residue_list):
            for j, resj in enumerate(chain._residue_list):
                if i == j: continue
                pos_i = resi.pos
                pos_j = resj.pos
                coord_i = np.array([pos_i[0], pos_i[1], pos_i[2]])
                coord_j = np.array([pos_j[0], pos_j[1], pos_j[2]])
                dist = np.linalg.norm(coord_i - coord_j)
                                    
                if dist <= cutoff:
                    # check for bonds
                    is_covalent = 1.0 if resi.covalent_bond(resj) else 0.0
                    is_hbond = 1.0 if resi.hydrogen_bond(
                        resj, 
                        hydrogen_present=False
                    ) else 0.0
                    edge_index.append([i, j])
                    edge_attr.append([dist, 
                                      is_covalent, 
                                      is_hbond]) 
                
        return np.array(edge_index).T, np.array(edge_attr)
    
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

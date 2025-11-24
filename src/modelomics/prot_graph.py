# created by clay 07/01/25
'''
create a graph representation of a protein via the modelomics.PDB object
'''
import numpy as np
from scipy.spatial import KDTree

from .utils.topology import vdwRadii

# properties of the atom to encode to numerical variables
element_to_id = {}
residue_to_id = {}

# a function to encode categorical variables into integers
def encode_categories(items, existing_dict=None):
    '''
    take a group of items (like atom names or residue names) and
    tokenize them
    '''
    # create a new dictionary if none provided
    if existing_dict is None:
        existing_dict = {}
    # a list of ID's
    ids = []
    for item in items:
        # is this item isnt already seen add a new encoding
        if item not in existing_dict:
            existing_dict[item] = len(existing_dict)
            # otherwise add it's encoding
        ids.append(existing_dict[item])
    # return the array of encoded variables and the updated dictionary 
    return np.array(ids), existing_dict

# function to turn a structure into features for the graph
def parse_pdb(pdb, ca_only=False, chains = False):
    """
    parse atoms from a PDB structure object.
    """
    positions = []
    numbers = []
    radii = []
    elements = []
    residues = []
    resnums = []

    if ca_only:
        if chains:
            for i, atom in enumerate(pdb.atoms):
                if atom.name == "CA" and atom.chain in chains:
                    positions.append(np.array([atom.x, atom.y, atom.z]))
                    numbers.append(atom.number)
                    radii.append(vdwRadii[atom.element])
                    elements.append(atom.element)
                    residues.append(atom.residue)
                    resnums.append(atom.resnum)
        else: 
            for i, atom in enumerate(pdb.atoms):
                if atom.name == "CA":
                    positions.append(np.array([atom.x, atom.y, atom.z]))
                    numbers.append(atom.number)
                    radii.append(vdwRadii[atom.element])
                    elements.append(atom.element)
                    residues.append(atom.residue)
                    resnums.append(atom.resnum)
    else:
        if chains:
            for i, atom in enumerate(pdb.atoms):
                if atom.chain in chains:
                    positions.append(np.array([atom.x, atom.y, atom.z]))
                    numbers.append(atom.number)
                    radii.append(vdwRadii[atom.element])
                    elements.append(atom.element)
                    residues.append(atom.residue)
                    resnums.append(atom.resnum)
        else: 
            for i, atom in enumerate(pdb.atoms):
                positions.append(np.array([atom.x, atom.y, atom.z]))
                numbers.append(atom.number)
                radii.append(vdwRadii[atom.element])
                elements.append(atom.element)
                residues.append(atom.residue)
                resnums.append(atom.resnum)

    # return the features 
    return (
        np.array(positions),
        np.array(numbers),
        np.array(radii),
        elements,
        residues,
        np.array(resnums)
    )

# a function to construct the edges of the graph using a kdtree for efficent
# neighbor searching
def build_edges(positions, cutoff=15.0):
    ''' 
    build the kdtree from the locations of the atoms'''
    # construct the kdtree based on the coordinates of the atoms
    tree = KDTree(positions)
    # get neighbors from some cutoff radius (all pairs with distance < r)
    pairs = tree.query_pairs(r=cutoff, output_type='ndarray')
    # turn the pairlist into forward and backwar pairs and transpose for graph
    edges = np.vstack([pairs, pairs[:, [1, 0]]]).T
    return edges

def pdb_to_pyg(pdb, ca_only=False, chains = None, cutoff=15.0):
    '''
    create a torch geometric graph representation from a modelomics.PDB object
    '''
    import torch
    from torch_geometric.data import Data

    pos_np, numbers, radii, elements, residues, resnums = parse_pdb(pdb, 
                                                                    ca_only, 
                                                                    chains)

    # encode categorical features to integers
    global element_to_id, residue_to_id
    element_ids, element_to_id = encode_categories(elements, element_to_id)
    residue_ids, residue_to_id = encode_categories(residues, residue_to_id)

    # turn the features into a matrix for torch graph
    features_np = np.stack([
        numbers,
        radii,
        element_ids,
        residue_ids, 
        resnums
    ], axis=1)

    # construct the feature matrix and position matrix
    x = torch.tensor(features_np, dtype=torch.float)
    pos = torch.tensor(pos_np, dtype=torch.float)

    # calculate the edges from the pairs that are close together with some r
    edge_index_np = build_edges(pos_np, cutoff=cutoff)
    # node indices must be integers
    edge_index = torch.tensor(edge_index_np, dtype=torch.long)

    # construct the graph
    return Data(x=x, edge_index=edge_index, pos=pos)
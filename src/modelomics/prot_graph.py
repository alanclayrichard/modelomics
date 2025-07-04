import numpy as np
from scipy.spatial import KDTree
import os


# Atomic radii lookup (values from a ucsf chimera page)
AtomicRadii = {
    "H": 1.0, "C": 1.7, "N": 1.625, "O": 1.48, "S": 1.782
}

# properties of the atom to encode to numerical variables
element_to_id = {}
residue_to_id = {}

# function to turn a structure into features for the graph
def parse_structure(structure, chain_id=None):
    """Parse atoms from a Biopython structure object."""
    positions = []
    serials = []
    radii = []
    elements = []
    residues = []

    # iterate through the Biopython models
    for model in structure:
        # iterate through the chains in the model
        for chain in model:
            if chain_id and chain.id != chain_id:
                continue
            for residue in chain:
                resname = residue.get_resname()
                for atom in residue:
                    # get the coordinates of this atom
                    coord = atom.coord
                    # the serial number from this atom
                    serial = atom.serial_number
                    # the element from this atom's name
                    element = atom.element.strip().capitalize()
                    # append the data to the lists
                    # Store features
                    positions.append(coord)
                    serials.append(serial)
                    radii.append(AtomicRadii.get(element, 0.0))
                    elements.append(element)
                    residues.append(resname)

    # return the features 
    return (
        np.array(positions),
        np.array(serials),
        np.array(radii),
        elements,
        residues
    )

# a function to load a pdb file into a biopython structure then parse
def parse_pdb(pdb_file, chain_id=None):
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pdb_struct", pdb_file)
    return parse_structure(structure, chain_id)

# a function to load a cif file into a biopython structure then parse
def parse_cif(cif_file, chain_id=None):
    from Bio.PDB import MMCIFParser
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("cif_struct", cif_file)
    return parse_structure(structure, chain_id)

# a function to encode categorical variables into integers
def encode_categories(items, existing_dict=None):
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

# a function to construct the edges of the graph using a kdtree for efficent
# neighbor searching
def build_edges(positions, cutoff=15.0):
    # construct the kdtree based on the coordinates of the atoms
    tree = KDTree(positions)
    # get neighbors from some cutoff radius (all pairs with distance < r)
    pairs = tree.query_pairs(r=cutoff, output_type='ndarray')
    # turn the pairlist into forward and backwar pairs and transpose for graph
    edges = np.vstack([pairs, pairs[:, [1, 0]]]).T
    return edges

def structure_to_pyg(file_path, chain=None, cutoff=15.0):
    import torch
    from torch_geometric.data import Data
    ext = os.path.splitext(file_path)[-1].lower()

    # parse the file according to its type (pdb and cif supported)
    if ext == ".pdb":
        pos_np, serials_np, radii_np, elements, residues = parse_pdb(
            file_path, 
            chain_id=chain
        )
    elif ext == ".cif":
        pos_np, serials_np, radii_np, elements, residues = parse_cif(
            file_path, 
            chain_id=chain
        )
    else:
        raise ValueError("Unsupported file format: must be .pdb or .cif")


    # encode categorical features to integers
    global element_to_id, residue_to_id
    element_ids, element_to_id = encode_categories(elements, element_to_id)
    residue_ids, residue_to_id = encode_categories(residues, residue_to_id)

    # turn the features into a matrix for torch graph (4 features per atom)
    features_np = np.stack([
        serials_np,
        radii_np,
        element_ids,
        residue_ids
    ], axis=1)  # shape: [num_atoms, 4]

    # construct the feature matrix and position matrix
    x = torch.tensor(features_np, dtype=torch.float)
    pos = torch.tensor(pos_np, dtype=torch.float)

    # calculate the edges from the pairs that are close together with some r
    edge_index_np = build_edges(pos_np, cutoff=cutoff)
    # node indices must be integers
    edge_index = torch.tensor(edge_index_np, dtype=torch.long)

    # construct the graph
    return Data(x=x, edge_index=edge_index, pos=pos)

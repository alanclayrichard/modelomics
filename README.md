# modelomics

`modelomics` is a comprehensive Python package for structural biology and protein design, providing state-of-the-art computational tools for protein analysis, graph neural networks, molecular surface analysis, and machine learning applications. Built with PyTorch Geometric integration for modern deep learning on protein structures.

## Current Features

### Hierarchical Protein Structure Representation
- **PDB parsing**: Complete PDB file parser with support for ATOM and HETATM records
- **Hierarchical organization**: Four-level structure hierarchy (PDB → Chain → Residue → Atom)
- **Flexible indexing**: Access structures by name or integer index at any level
- **Structure cleaning**: Remove waters, heteroatoms, and filter by specific chains

### Atom-Level Processing  
- **Comprehensive atom properties**: Coordinates, van der Waals radii, element types, occupancy, B-factors
- **Insertion code support**: Full support for PDB insertion codes in residue numbering

### Chain-Level Operations
- **Dictionary and list access**: Access residues by `("name", number, "insertion")` tuples or integer index

### Advanced Residue-Level Analysis
- **Comprehensive atom organization**: Residues contain atoms as both dictionary (by name) and list (by index)
- **Structural analysis methods**: Alpha carbon positioning, geometric calculations, coordinate extraction
- **Interaction detection**: Covalent bond detection using atomic radii and distance criteria
- **Hydrogen bond analysis**: Geometric hydrogen bond detection with donor-acceptor distance and angle constraints

### Protein Sequence Analysis and Embeddings
- **ESM-C language model embeddings**: Generate sequence embeddings from pretrained models
- **Masked sequence analysis**: Create embeddings with specific residue positions masked for mutation effect studies
- **Embedding integration**: Integration of sequence embeddings into graph neural network node features

### Advanced Graph Neural Network Representations
- **PyTorch Geometric integration**: Native `torch_geometric.data.Data` objects for GNN training
- **Rich residue-level graphs**: Features combining sequence embeddings, physicochemical properties, SASA, and normal mode analysis
- **Multi-scale edge features**: Distance, covalent bonds, hydrogen bonds, and hydrophobic contacts
- **ESM-C sequence embeddings**: 960-dimensional protein language model representations
- **Structural features**: SASA analysis (5D), normal mode fluctuations, physicochemical properties (5D)
- **Efficient bond detection**: Fast algorithms for covalent bonding, hydrogen bonding

### Molecular Surface Analysis and SASA Calculations
- **Fibonacci sphere sampling**: Mathematically optimal point distribution on atomic surfaces (92 points per atom)
- **Van der Waals surface modeling**: Accurate molecular surface generation with customizable probe radius (1.4Å default)
- **SASA calculations**: Solvent-accessible surface area calculations with atomic-level precision
- **Buried point removal**: KDTree-based neighbor search algorithms for efficient buried point identification
- **Typed SASA analysis**: Surface area breakdown by polarity (polar/nonpolar/positive/negative/total)
- **Residue-level SASA**: Aggregate surface area calculations at the residue level for binding site analysis
- **Visualization export**: Export point clouds in XYZ format with atomic coloring for molecular visualization

## Installation

### Prerequisites

- Python 3.10 or higher

### Install with pip

```bash
pip install git+https://github.com/alanclayrichard/modelomics.git
```

## Example Usage

### Hierarchical Structure Access

```python
from modelomics.PDB import PDB

# Load a PDB structure
pdb = PDB("protein.pdb")

# Access by chain name or index
chain_A = pdb['A']          # Access by name
first_chain = pdb[0]        # Access by index

# Access residues within chains
residue = chain_A['GLN', 1, '']    # Access by (name, number, insertion)
residue = chain_A[0]               # Access by index

# Access atoms within residues  
atom = residue['CA']        # Access by atom name
atom = residue[0]          # Access by index

# Complete hierarchical access
ca_atom = pdb['A']['GLN', 1, '']['CA']  # Chain A, GLN 1, CA atom
```

### Structure Processing and Filtering

```python
from modelomics.PDB import PDB

# Load and clean structure
pdb = PDB("protein.pdb")
print(f"Original: {len(pdb.atoms)} atoms in {len(pdb.chains)} chains")

# Remove waters and heteroatoms
pdb.strip_waters()  # Remove HOH residues
pdb.strip_hetatm()  # Keep only ATOM records
print(f"Cleaned: {len(pdb.atoms)} atoms")

# Access cleaned structure
for chain_id, chain in pdb.chains.items():
    print(f"Chain {chain_id}: {len(chain.residues)} residues")
    for residue in chain._residue_list[:3]:  # First 3 residues
        print(f"  {residue.name}{residue.number}: {len(residue.atoms)} atoms")
```

### Sequence Analysis and Embeddings
```python
from modelomics import sequences

# Extract sequences from structures
seq_cif = sequences.sequence_from_cif("protein.cif", chain='B')
seq_pdb = sequences.sequence_from_pdb("protein.pdb", chain='A')
print("Sequence:", seq_pdb)

# Generate ESM-C protein language model embeddings
embedding = sequences.embed_sequence(seq_pdb)
print("Embedding shape:", embedding.shape)  # [1, seq_len+2, hidden_dim]

# Create masked embeddings for mutation analysis
masked_embedding = sequences.embed_sequence_with_mask(seq_pdb, mask_positions=[5, 10])
print("Masked embedding shape:", masked_embedding.shape)
```

### Advanced Graph Neural Network Representations

```python
from modelomics import PDB, Graph
import torch

# Load and prepare protein structure
pdb = PDB("protein.pdb")
pdb.strip_waters()
pdb.rename_to_charmm()  # Standardize atom naming

# Create rich residue-level graph with comprehensive features
graph = Graph(pdb["A"], cutoff=15.0)
print(f"Graph: {graph.num_nodes} residues, {graph.num_edges} edges")
print(f"Node features: {graph.x.shape} (971D: 960 ESM + 5 physicochemical + 5 SASA + 1 NMA)")
print(f"Edge features: {graph.edge_attr.shape} (distance + covalent + H-bonds + hydrophobic)")

# Access rich structural features
print("\nNode feature breakdown:")
print(f"- ESM-C embeddings: {graph.x[:, :960].shape}")
print(f"- Physicochemical properties: {graph.x[:, 960:965].shape}")
print(f"- SASA features: {graph.x[:, 965:970].shape}")
print(f"- Normal mode analysis: {graph.x[:, 970:971].shape}")

# Examine edge features for molecular interactions
print(f"\nEdge features: [distance, covalent_bond, hydrogen_bond, hydrophobic_contact]")
print(f"Sample edges: {graph.edge_attr[:5]}")

# Ready for PyTorch Geometric neural networks
print(f"\nPyTorch Geometric Data object: {type(graph)}")
print(f"Device ready: {graph.x.device}")
```

### Molecular Surface Analysis and SASA Calculations

```python
from modelomics import PDB, VDWPointCloud

# Load and prepare structure
pdb = PDB("protein.pdb")
pdb.strip_waters()
chain_atoms = pdb['A'].atoms

# Generate high-quality molecular surface
pc = VDWPointCloud(chain_atoms, N=92, probe_radius=1.4)
pc.populate()  # Generate Fibonacci spheres on all atoms
pc.remove_buried()  # KDTree-based buried point removal

# Comprehensive SASA analysis
total_sasa = pc.sasa_by_type()
print("Global SASA breakdown:")
for sasa_type, value in total_sasa.items():
    print(f"  {sasa_type}: {value:.2f} Ų")

# Residue-specific SASA for binding site analysis
binding_site_residues = [pdb['A'][i] for i in range(10, 20)]  # Example binding site
site_sasa = pc.residue_sasa_by_type(binding_site_residues)
print(f"\nBinding site SASA: {site_sasa['total']:.2f} Ų")
print(f"Hydrophobic surface: {site_sasa['nonpolar']:.2f} Ų")

# Per-atom SASA calculations
for i, atom in enumerate(chain_atoms[:5]):
    atom_sasa = pc.get_atom_sasa(i)
    print(f"{atom.residue}{atom.resnum} {atom.name}: {atom_sasa:.2f} Ų")

# Export surface for molecular visualization
pc.write_points("surface.xyz")
total_points = sum(len(pts) for pts in pc.points)
print(f"\nExported {total_points} surface points to surface.xyz")
```

### Residue-Level Structural Analysis

```python
from modelomics import PDB

# Load and prepare protein complex
pdb = PDB("complex.pdb")
pdb.strip_waters()
pdb.rename_to_charmm()

# Advanced residue interaction analysis
chain_a = pdb['A']
chain_b = pdb['B']

# Find interface residues with molecular interactions
interface_pairs = []
for res_a in chain_a._residue_list:
    for res_b in chain_b._residue_list:
        # Check all types of molecular interactions
        if res_a.covalent_bond(res_b):
            interface_pairs.append((res_a, res_b, "covalent"))
        elif res_a.hydrogen_bond(res_b):
            interface_pairs.append((res_a, res_b, "hydrogen_bond"))

print(f"Found {len(interface_pairs)} interface interactions:")
for res1, res2, interaction_type in interface_pairs[:10]:  # First 10
    pos1 = res1.pos
    pos2 = res2.pos
    distance = ((pos1[0]-pos2[0])**2 + (pos1[1]-pos2[1])**2 + (pos1[2]-pos2[2])**2)**0.5
    print(f"{res1.name}{res1.number} - {res2.name}{res2.number}: {interaction_type} ({distance:.2f}Å)")

```


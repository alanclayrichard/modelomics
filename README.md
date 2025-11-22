# modelomics

`modelomics` is a Python package for structural biology and protein design, providing computational tools for protein analysis, visualization, and machine learning applications.

## Current Features

### Hierarchical Protein Structure Representation
- **PDB parsing**: Complete PDB file parser with support for ATOM and HETATM records
- **Hierarchical organization**: Four-level structure hierarchy (PDB → Chain → Residue → Atom)
- **Flexible indexing**: Access structures by name or integer index at any level
- **Structure cleaning**: Remove waters, heteroatoms, and filter by specific chains

### Atom-Level Processing  
- **Comprehensive atom properties**: Coordinates, van der Waals radii, element types, occupancy, B-factors
- **Insertion code support**: Full support for PDB insertion codes in residue numbering
- **Polarity classification**: Built-in methods to identify polar/nonpolar atoms
- **Element validation**: Automatic element detection and van der Waals radius assignment

### Chain-Level Operations
- **Automatic chain detection**: Efficient grouping of atoms by chain identifier
- **Residue organization**: Smart grouping of atoms into residues by number, name, and insertion code
- **Dictionary and list access**: Access residues by `("name", number, "insertion")` tuples or integer index
- **Chain filtering**: Maintain chain integrity during structure manipulation

### Residue-Level Management
- **Atom organization**: Residues contain atoms as both dictionary (by name) and list (by index)
- **Flexible atom access**: Get atoms by name (`residue["CA"]`) or index (`residue[0]`)
- **Insertion code handling**: Proper support for PDB insertion codes in residue identification
- **Consistency validation**: Automatic verification of residue-atom relationships

### Sequence Analysis  
- **Multi-format sequence extraction**: Extract amino acid sequences from PDB and CIF files
- **Protein language model embeddings**: Generate sequence embeddings using Meta's ESM-C models
- **Masked sequence embeddings**: Create embeddings with specific residue positions masked for mutation studies
- **Lazy model loading**: Efficient memory management with on-demand model loading

### Graph Neural Network Integration
- **PyTorch Geometric compatibility**: Convert protein structures to graph neural network format
- **Multi-level graphs**: Support for atom-level or residue-level (CA-only) graph representations
- **Spatial connectivity**: Automatic edge construction based on distance cutoffs using KDTree algorithms
- **Rich node features**: Atomic number, van der Waals radius, element ID, residue type, residue number
- **Flexible filtering**: Filter by specific chains or atom types during graph construction

### 3D Point Cloud Generation
- **Fibonacci sphere sampling**: Mathematically optimal point distribution on atomic surfaces
- **Van der Waals surface modeling**: Generate molecular surfaces with probe radius support
- **Buried point removal**: KDTree-based algorithms to identify solvent-accessible surface points
- **Efficient sphere generation**: Pre-computed unit vectors for fast point cloud generation
- **Visualization export**: Export point clouds in XYZ format for molecular visualization tools

### Utilities and Data
- **Comprehensive van der Waals radii**: Database covering 100+ elements from standard references
- **CPK color schemes**: Standard atomic color schemes for molecular visualization
- **Categorical encoding**: Efficient encoding of molecular features for machine learning applications
- **Memory optimization**: Efficient data structures and algorithms for large protein systems

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

### Graph Neural Network Representations

```python
from modelomics.PDB import PDB
from modelomics.prot_graph import pdb_to_pyg

# Create atom-level graph
pdb = PDB("protein.pdb")
pdb.strip_waters()

# Full atom graph
atom_graph = pdb_to_pyg(pdb, ca_only=False, chains=['A'], cutoff=15.0)
print(f"Atom graph: {atom_graph.x.shape[0]} nodes, {atom_graph.edge_index.shape[1]} edges")

# CA-only residue graph  
ca_graph = pdb_to_pyg(pdb, ca_only=True, chains=['A'], cutoff=15.0)
print(f"CA graph: {ca_graph.x.shape[0]} nodes, {ca_graph.edge_index.shape[1]} edges")

# Access node features: [atomic_number, vdw_radius, element_id, residue_id, residue_number]
print("Node features:", atom_graph.x[:3])  # First 3 nodes
print("3D positions:", atom_graph.pos[:3])  # First 3 positions
```

### Van der Waals Surface Point Clouds

```python
from modelomics.PDB import PDB
from modelomics.vdw_point_cloud import vdwPointCloud

# Load and prepare structure
pdb = PDB("protein.pdb")
pdb.strip_waters()
pdb.strip_hetatm()

# Generate molecular surface point cloud
point_cloud = vdwPointCloud(pdb.atoms, N=92, probe_radius=1.4)
print(f"Point cloud for {len(point_cloud.atoms)} atoms")

# Generate surface points
point_cloud.populate()
print(f"Generated {sum(len(pts) for pts in point_cloud.points)} surface points")

# Remove buried points to get solvent-accessible surface
point_cloud.remove_buried()
accessible_points = sum(len(pts) for pts in point_cloud.points)
print(f"Solvent-accessible points: {accessible_points}")

# Export for visualization in Ovito or other tools
point_cloud.write_points("protein_surface.xyz")
```

### Advanced Usage Patterns

```python
from modelomics.PDB import PDB

# Load structure
pdb = PDB("complex.pdb")

# Analyze multi-chain complex
for chain_id in pdb.chains:
    chain = pdb[chain_id]
    print(f"Chain {chain_id}: {len(chain._residue_list)} residues")
    
    # Find interface residues (example with distance cutoff)
    if len(pdb.chains) > 1:
        other_chains = [c for c in pdb.chains if c != chain_id]
        # Analysis logic here...

# Work with specific residue types
for residue in pdb['A']._residue_list:
    if residue.name in ['ARG', 'LYS', 'HIS']:  # Positive residues
        print(f"Positive residue: {residue}")
        
        # Access specific atoms
        if 'CA' in residue.atoms:
            ca = residue['CA']
            print(f"  CA coordinates: ({ca.x:.2f}, {ca.y:.2f}, {ca.z:.2f})")
```


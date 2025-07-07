# modelomics

`modelomics` is a Python package for working with various models and tools for protein design.

## Current Features

- Easy integration with Biopython's PDB parsing utilities.
- Extract sequences from CIF structure files.
- Generate sequence embeddings using ESM models.
- Generate graph representations of protein structures using PyTorch Geometric.


## Installation

### Prerequisites

- Python 3.10 or higher

### Install with pip

```bash
pip install git+https://github.com/alanclayrichard/modelomics.git
```

## Example Usage

This package provides tools for building protein graphs from PDB/CIF files and embedding sequences using pretrained models.

### Build Atom- or Residue-Level Graphs

```python
from modelomics.atom_graph import AtomGraphBuilder
from modelomics.residue_graph import ResidueGraphBuilder

# Build atom-level graph from a PDB file
atom_builder = AtomGraphBuilder()
atom_graph = atom_builder.build("pdbs/3ux9.pdb")

# Build residue-level graph from the same structure
residue_builder = ResidueGraphBuilder()
residue_graph = residue_builder.build("pdbs/3ux9.cif")

# Access node and edge info
print("Atom graph nodes:", atom_graph.x.shape)
print("Residue graph edges:", residue_graph.edge_index.shape)
```
### Extract and Embed Sequences
```python

from modelomics import sequences

# Extract amino acid sequence from a chain
seq = sequences.sequence_from_cif("pdbs/3ux9.cif", chain='B')
print("Sequence:", seq)

# Embed the sequence using a pretrained language model (e.g., ESM)
embedding = sequences.embed_sequence(seq)
print("Embedding shape:", embedding.shape)  # [1, seq_len+2, hidden_dim]

# Optionally embed with site masking (e.g., mutate site 1 and 2)
masked_embedding = sequences.embed_sequence_with_mask(seq, mask_sites=[1, 2])
```


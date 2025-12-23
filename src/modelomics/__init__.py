"""
modelomics: A Python package for structural biology and protein design
"""

__author__ = "A. Clay Richard"
__email__ = "alanclayrichard@gmail.com"

# Lazy imports to speed up package loading
def __getattr__(name):
    if name == "PDB":
        from . import pdb
        return pdb.PDB
    elif name == "CIF":
        from . import cif
        return cif.CIF
    elif name == "Atom":
        from . import atom
        return atom.Atom
    elif name == "Chain":
        from . import chain
        return chain.Chain
    elif name == "Residue":
        from . import residue
        return residue.Residue
    elif name == "VDWPointCloud":
        from . import vdw_point_cloud
        return vdw_point_cloud.VDWPointCloud
    elif name == "sequence_embeddings":
        from . import sequence_embeddings
        return sequence_embeddings
    elif name == "Graph":
        from . import graph
        return graph.Graph
    else:
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")

__all__ = [
    "PDB",
    "CIF",
    "Atom", 
    "Chain",
    "Residue",
    "VDWPointCloud", 
    "sequence_embeddings",
    "Graph",
    "__version__",
    "__author__",
    "__email__"
]

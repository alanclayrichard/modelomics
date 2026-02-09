# created by clay 02/09/26
'''
CHARMM atom names for the 20 standard amino acids. these follow the naming
conventions applied by atom.rename_to_charmm (leading digit moved to end,
backbone H -> HN, ILE CD1 -> CD, terminal O -> OT1/OT2, etc.)
'''

import numpy as np

# sorted unique CHARMM atom names across all 20 standard amino acids
# backbone + sidechain + terminal atoms
charmm_atom_names = [
    # backbone
    "C", "CA", "CB",
    # sidechain carbons
    "CD", "CD1", "CD2",
    "CE", "CE1", "CE2", "CE3",
    "CG", "CG1", "CG2",
    "CH2", "CZ", "CZ2", "CZ3",
    # alpha hydrogens
    "HA", "HA1", "HA2",
    # beta hydrogens
    "HB", "HB1", "HB2", "HB3",
    # delta hydrogens
    "HD1", "HD2", "HD3",
    "HD11", "HD12", "HD13",
    "HD21", "HD22", "HD23",
    # epsilon hydrogens
    "HE", "HE1", "HE2", "HE3",
    "HE21", "HE22",
    # gamma hydrogens
    "HG", "HG1", "HG2",
    "HG11", "HG12", "HG13",
    "HG21", "HG22", "HG23",
    # ring / guanidinium hydrogens
    "HH", "HH2",
    "HH11", "HH12", "HH21", "HH22",
    # backbone amide hydrogen (CHARMM: HN not H)
    "HN",
    # n-terminal hydrogens
    "HT1", "HT2", "HT3",
    # zeta hydrogens
    "HZ", "HZ1", "HZ2", "HZ3",
    # backbone nitrogen
    "N",
    # sidechain nitrogens
    "ND1", "ND2", "NE", "NE1", "NE2",
    "NH1", "NH2", "NZ",
    # backbone oxygen
    "O",
    # sidechain oxygens
    "OD1", "OD2", "OE1", "OE2",
    "OG", "OG1", "OH",
    # c-terminal oxygens
    "OT1", "OT2",
    # sulfur
    "SD", "SG",
]

charmm_atom_name_to_idx = {name: i for i, name in enumerate(charmm_atom_names)}
CHARMM_ATOM_ONEHOT = np.eye(len(charmm_atom_names))

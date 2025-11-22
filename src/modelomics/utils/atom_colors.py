# created by clay 11/22/25
'''
atom rgb colors for use in various visualization
'''
import random

from vdw_radii import vdwRadii

CPK_colors = {
    "H":  (1.0, 1.0, 1.0),      # white
    "C":  (0.5, 0.5, 0.5),      # dark gray
    "N":  (0.0, 0.0, 1.0),      # blue
    "O":  (1.0, 0.0, 0.0),      # red
    "F":  (0.0, 1.0, 0.0),      # green
    "CL": (0.0, 1.0, 0.0),      # green
    "BR": (0.6, 0.13, 0.0),     # dark red/brown
    "I":  (0.4, 0.0, 0.6),      # purple
    "HE": (0.85, 1.0, 1.0),     # cyan
    "NE": (0.7, 0.89, 0.96),    # light blue
    "AR": (0.5, 0.82, 0.89),    # light blue
    "XE": (0.25, 0.61, 0.81),   # medium blue
    "KR": (0.36, 0.72, 0.82),   # cyan-ish
    "P":  (1.0, 0.5, 0.0),      # orange
    "S":  (1.0, 1.0, 0.0),      # yellow
    "B":  (1.0, 0.7, 0.7),      # pink
    "LI": (0.8, 0.5, 1.0),      # violet
    "NA": (0.67, 0.36, 0.95),   # purple
    "K":  (0.56, 0.25, 0.83),   # dark purple
    "CA": (0.24, 1.0, 0.0),     # green
    "MG": (0.12, 1.0, 0.0),     # green
    "FE": (1.0, 0.6, 0.0),      # orange
    "CU": (0.78, 0.5, 0.2),     # brown
    "ZN": (0.49, 0.5, 0.69),    # purple-gray
}

atom_colors = {}
for element in vdwRadii:
    if element in CPK_colors:
        atom_colors[element] = CPK_colors[element]
    else:
        # Assign a random pastel color if unknown
        atom_colors[element] = tuple(random.uniform(0.5, 1.0) for _ in range(3))

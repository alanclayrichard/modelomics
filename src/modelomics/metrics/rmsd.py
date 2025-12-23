# created by clay
'''
calculate rmsd of sets of atoms
'''

import numpy as np

def rmsd_atoms(atoms_1, atoms_2):
    '''
    get rmsd of two sets of atoms assuming theyre in the 
    right order
    '''
    coords_1 = np.array([[atom.x, atom.y, atom.z] for atom in atoms_1])
    coords_2 = np.array([[atom.x, atom.y, atom.z] for atom in atoms_2])
    return rmsd(coords_1, coords_2)

def rmsd(x1, x2):
    diff = x1 - x2
    return np.sqrt((diff * diff).sum(axis=1).mean())
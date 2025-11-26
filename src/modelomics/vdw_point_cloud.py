# created by clay 07/01/25
"""
efficient generation of Fibonacci spheres for atoms,
with buried point removal using KDTree neighbor search.
"""

import numpy as np
from scipy.spatial import KDTree

from .utils.atom_colors import atom_colors

class VDWPointCloud:
    def __init__(self, atoms, N=92, probe_radius=1.4):
        """
        initialize a point cloud from a list of modelomics.atom objects
        """
        self.atoms = atoms
        self.N = N
        self.probe_radius = probe_radius

        # store points as a list of arrays (variable-length)
        self.points = []

        # precompute Fibonacci sphere directions (unit vectors)
        golden_ratio = (1 + np.sqrt(5)) / 2
        i = np.arange(N)
        polar = np.arccos(1 - 2*(i + 0.5)/N)
        azimuth = 2 * np.pi * i / golden_ratio

        self.unit_dirs = np.column_stack([
            np.sin(polar) * np.cos(azimuth),
            np.sin(polar) * np.sin(azimuth),
            np.cos(polar)
        ])

    def __str__(self):
        return (f"Point cloud with {len(self.atoms)} atoms, "
               f"spheres of {self.N} points each")

    def generate_sphere(self, atom, probe_radius=0):
        """
        Generate the Fibonacci sphere for a given atom
        """
        R_eff = atom.radius + probe_radius
        center = np.array([atom.x, atom.y, atom.z])
        return center + R_eff * self.unit_dirs

    def populate(self):
        """
        populate the point cloud for all atoms
        """
        self.points = [self.generate_sphere(atom, self.probe_radius) 
                       for atom in self.atoms]

    def remove_buried(self):
        """
        remove points buried by neighboring atoms using KDTree
        """
        atom_positions = np.array([[atom.x, atom.y, atom.z] 
                                   for atom in self.atoms])
        atom_radii = np.array([atom.radius + self.probe_radius 
                               for atom in self.atoms])

        # build KDTree for neighbor search
        tree = KDTree(atom_positions)

        for i, atom in enumerate(self.atoms):
            points = self.points[i]  # shape (N,3)

            # find neighboring atoms within max radius distance
            neighbors = tree.query_ball_point([atom.x, 
                                               atom.y, 
                                               atom.z], 
                                               r=np.max(atom_radii)*2)
            neighbors = [idx for idx in neighbors if idx != i]

            neighbor_positions = atom_positions[neighbors]
            neighbor_radii = atom_radii[neighbors]

            # vectorized distance check
            diff = points[:, np.newaxis, :]-neighbor_positions[np.newaxis, :, :]
            dist_sq = np.sum(diff**2, axis=2)
            buried = np.any(dist_sq < neighbor_radii[np.newaxis, :]**2, axis=1)

            # keep only non-buried points
            self.points[i] = points[~buried]

    def get_atom_sasa(self, atom_idx):
        """
        calculate SASA for a specific atom based on exposed points
        """
        if atom_idx >= len(self.points):
            return 0.0
        
        atom = self.atoms[atom_idx]
        exposed_points = len(self.points[atom_idx])
        total_points = self.N
        
        # calculate surface area of the atom's vdw sphere
        R_eff = atom.radius + self.probe_radius
        total_surface_area = 4 * np.pi * R_eff**2
        
        # fraction exposed * total surface area
        return (exposed_points / total_points) * total_surface_area

    def get_residue_sasa(self, residue):
        """
        calculate total SASA for a residue by summing its atoms
        """
        total_sasa = 0.0
        for atom in residue._atom_list:
            # find this atom in our atom list
            for i, cloud_atom in enumerate(self.atoms):
                if (atom.x == cloud_atom.x and 
                    atom.y == cloud_atom.y and 
                    atom.z == cloud_atom.z):
                    total_sasa += self.get_atom_sasa(i)
                    break
        return total_sasa

    def sasa_by_type(self):
        """
        calculate total and typed SASA for entire protein based on atom polarity
        """
        polar_sasa = 0.0
        nonpolar_sasa = 0.0
        positive_sasa = 0.0
        negative_sasa = 0.0
        total_sasa = 0.0

        for i, atom in enumerate(self.atoms):
            atom_sasa = self.get_atom_sasa(i)
            total_sasa += atom_sasa

            # classify by atom polarity (from atom.py)
            if atom.polar:
                polar_sasa += atom_sasa
            else:
                nonpolar_sasa += atom_sasa
                
            # classify by formal charge
            if hasattr(atom, 'charge'):
                if atom.charge > 0:
                    positive_sasa += atom_sasa
                elif atom.charge < 0:
                    negative_sasa += atom_sasa

        return {
            'total': total_sasa,
            'polar': polar_sasa,
            'nonpolar': nonpolar_sasa,
            'positive': positive_sasa,
            'negative': negative_sasa
        }

    def residue_sasa_by_type(self, residues):
        """
        calculate SASA by atom type for specific binding site residues
        """
        polar_sasa = 0.0
        nonpolar_sasa = 0.0
        positive_sasa = 0.0
        negative_sasa = 0.0
        total_sasa = 0.0

        for res in residues:
            try:
                if hasattr(res, 'is_missing') and res.is_missing:
                    continue
                    
                # calculate SASA for each atom in the residue
                for atom in res._atom_list:
                    # find this atom in our atom list
                    for i, cloud_atom in enumerate(self.atoms):
                        if (atom.x == cloud_atom.x and 
                            atom.y == cloud_atom.y and 
                            atom.z == cloud_atom.z):
                            
                            atom_sasa = self.get_atom_sasa(i)
                            total_sasa += atom_sasa
                            
                            # classify by atom polarity
                            if atom.polar:
                                polar_sasa += atom_sasa
                            else:
                                nonpolar_sasa += atom_sasa
                                
                            # classify by formal charge
                            if hasattr(atom, 'charge'):
                                if atom.charge > 0:
                                    positive_sasa += atom_sasa
                                elif atom.charge < 0:
                                    negative_sasa += atom_sasa
                            break
                            
            except KeyError:
                continue

        return {
            'total': total_sasa,
            'polar': polar_sasa,
            'nonpolar': nonpolar_sasa,
            'positive': positive_sasa,
            'negative': negative_sasa
        }

    def write_points(self, xyz_file):
        """
        write points to an xyz file for visualization
        """
        total_points = sum(len(p) for p in self.points) - 1
        with open(xyz_file, 'w') as f:
            f.write(f"{total_points}\n")
            for i, atom_points in enumerate(self.points):
                element = self.atoms[i].element
                color = atom_colors[element]
                for pt in atom_points:
                    f.write(f"{pt[0]} {pt[1]} {pt[2]} 0.1 "
                            f"{color[0]} {color[1]} {color[2]}\n")

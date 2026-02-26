# created by clay 02/26/26
'''
the Trajectory class — a list of PDB frames backed by mdtraj
'''


class Trajectory:
    def __init__(self, source, top=None):
        '''
        initialize a Trajectory object

        Parameters
        ----------
        source : mdtraj.Trajectory, str, or list of PDB
            - mdtraj.Trajectory: wrap directly
            - str: load from file (DCD, XTC, TRR, PDB, ...)
            - list: build from a list of PDB objects
        top : str or None
            topology file, required when source is a trajectory file
        '''
        import mdtraj
        self._cache = {}

        if isinstance(source, mdtraj.Trajectory):
            self._traj = source
        elif isinstance(source, str):
            self._traj = mdtraj.load(source, top=top)
        elif isinstance(source, list):
            self._traj = self._join_pdbs(source)
        else:
            raise TypeError(f"cannot construct Trajectory from {type(source)}")

    def __str__(self):
        return (
            f"Trajectory with {self._traj.n_frames} frames, "
            f"{self._traj.n_atoms} atoms, "
            f"{self._traj.topology.n_residues} residues, "
            f"{self._traj.topology.n_chains} chains"
        )

    def __len__(self):
        return self._traj.n_frames

    def __getitem__(self, key):
        if isinstance(key, int):
            n = len(self)
            if key < 0:
                key += n
            if not 0 <= key < n:
                raise IndexError(f"frame index {key} out of range for trajectory with {n} frames")
            if key not in self._cache:
                self._cache[key] = self._frame_to_pdb(key)
            return self._cache[key]
        elif isinstance(key, slice):
            return Trajectory(self._traj[key])
        else:
            raise TypeError(f"invalid index type: {type(key)}")

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def n_frames(self):
        return self._traj.n_frames

    def n_atoms(self):
        return self._traj.n_atoms

    def n_residues(self):
        return self._traj.topology.n_residues

    def n_chains(self):
        return self._traj.topology.n_chains

    def xyz(self):
        '''coordinates in Angstroms, shape (n_frames, n_atoms, 3)'''
        return self._traj.xyz * 10.0

    def time(self):
        return self._traj.time

    def topology(self):
        return self._traj.topology

    def rmsd(self, reference=0, atom_indices=None):
        '''
        per-frame RMSD relative to a reference frame, in Angstroms

        Parameters
        ----------
        reference : int
            frame index to use as reference (default 0)
        atom_indices : array-like or None
            atom indices to include; None uses all atoms
        '''
        import mdtraj
        return mdtraj.rmsd(self._traj, self._traj, reference, atom_indices=atom_indices) * 10.0

    def superpose(self, reference=0, atom_indices=None):
        '''
        in-place superposition of all frames onto a reference frame

        Parameters
        ----------
        reference : int
            frame index to use as reference (default 0)
        atom_indices : array-like or None
            atom indices to align on; None uses all atoms
        '''
        self._traj.superpose(self._traj, reference, atom_indices=atom_indices)
        self._cache.clear()

    def save(self, filename, **kwargs):
        '''
        save the trajectory to a file

        Parameters
        ----------
        filename : str
            output path; format inferred from extension
        '''
        self._traj.save(filename, **kwargs)

    def _join_pdbs(self, pdb_list):
        import mdtraj
        import tempfile
        import os

        trajs = []
        tmp_files = []
        try:
            for pdb_obj in pdb_list:
                fd, path = tempfile.mkstemp(suffix='.pdb')
                os.close(fd)
                tmp_files.append(path)
                pdb_obj.write(path)
                trajs.append(mdtraj.load_pdb(path))
        finally:
            for path in tmp_files:
                try:
                    os.unlink(path)
                except OSError:
                    pass

        return mdtraj.join(trajs)

    def _frame_to_pdb(self, frame_idx):
        from . import pdb as p, atom as a

        top = self._traj.topology
        xyz = self._traj.xyz[frame_idx] * 10.0  # nm → Å

        lines = []
        for atom in top.atoms:
            x, y, z = xyz[atom.index]
            res = atom.residue
            chain = res.chain
            try:
                chain_id = chain.chain_id or chr(65 + chain.index % 26)
            except AttributeError:
                chain_id = chr(65 + chain.index % 26)
            element = atom.element.symbol if atom.element else ' '
            line = (
                f"{'ATOM':<6}"           # [0:6]   record type
                f"{atom.index+1:>5} "   # [6:11] + [11] space
                f"{atom.name:<4}"       # [12:16] atom name
                f" "                    # [16]    altloc
                f"{res.name:<3} "       # [17:20] resname + [20] space
                f"{chain_id:1}"         # [21]    chain ID
                f"{res.resSeq:>4}"      # [22:26] resSeq
                f" "                    # [26]    insertion code
                f"   "                  # [27:30] spaces
                f"{x:>8.3f}"            # [30:38] x
                f"{y:>8.3f}"            # [38:46] y
                f"{z:>8.3f}"            # [46:54] z
                f"  1.00"               # [54:60] occupancy
                f"  0.00"               # [60:66] temp_factor
                f"          "           # [66:76] spaces
                f"{element:>2}\n"       # [76:78] element
            )
            lines.append(line)

        pdb_obj = p.PDB.__new__(p.PDB)
        pdb_obj.filename = f"<mdtraj frame {frame_idx}>"
        pdb_obj.file_data = ""
        pdb_obj.missing_residues = {}
        pdb_obj.atoms = [a.Atom(line) for line in lines]
        pdb_obj.chains = {}
        pdb_obj._chain_list = []
        pdb_obj._build_chains()
        return pdb_obj

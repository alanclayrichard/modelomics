class vdwPointCloud:
    def __init__(self, atom_graph):
        self.atom_graph = atom_graph

    def __str__(self):
        return f"graph has {self.atom_graph.num_nodes} nodes"
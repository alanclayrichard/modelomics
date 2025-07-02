import unittest
import os
import torch
from torch_geometric.data import Data
from modelomics import prot_graph 

class TestGraphFromPDB(unittest.TestCase):
    # setup that runs once before all tests in this class
    @classmethod
    def setUpClass(cls):
        # get the location of the test file (could be pdb or cif)
        cls.test_dir = os.path.dirname(__file__)
        cls.protein_path = os.path.join(cls.test_dir, "pdbs", "3ux9.pdb")
        # build graph from file
        cls.data = prot_graph.structure_to_pyg(cls.protein_path, chain=None) 

    def test_graph_structure_from_pdb(self):
        data = self.data

        # check types and if theres something there
        self.assertIsInstance(data, Data)
        self.assertTrue(hasattr(data, 'x'))
        self.assertTrue(hasattr(data, 'edge_index'))
        self.assertTrue(hasattr(data, 'pos'))

        # check dimensions
        self.assertEqual(data.x.ndim, 2)
        self.assertEqual(data.x.shape[1], 4)
        self.assertEqual(data.pos.ndim, 2)
        self.assertEqual(data.edge_index.shape[0], 2)
        self.assertEqual(data.x.shape[0], data.pos.shape[0])

        # ensure its not empty
        self.assertGreater(data.x.shape[0], 0)
        self.assertGreater(data.edge_index.shape[1], 0)

    def test_no_isolated_nodes(self):
        data = self.data
        degrees = torch.zeros(data.num_nodes, dtype=torch.long)
        degrees.scatter_add_(0, 
                             data.edge_index[0], 
                             torch.ones(data.edge_index.size(1), 
                             dtype=torch.long))
        self.assertTrue(torch.all(degrees > 0), 
                        "All nodes should have at least one edge")

if __name__ == "__main__":
    unittest.main()

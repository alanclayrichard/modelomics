import unittest
import os
from modelomics.atom_graph import AtomGraphBuilder
from modelomics.residue_graph import ResidueGraphBuilder
from torch_geometric.data import Data

class TestAtomGraphBuilder(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.test_dir = os.path.dirname(__file__)
        cls.pdb_path = os.path.join(cls.test_dir, "pdbs", "3ux9.pdb")
        cls.cif_path = os.path.join(cls.test_dir, "pdbs", "3ux9.cif")
        cls.builder = AtomGraphBuilder()

    def test_build_graph_from_pdb(self):
        data = self.builder.build(filename=self.pdb_path)
        self._validate_graph(data)

    def test_build_graph_from_cif(self):
        data = self.builder.build(filename=self.cif_path)
        self._validate_graph(data)

    def _validate_graph(self, data):
        self.assertIsInstance(data, Data)
        self.assertTrue(hasattr(data, 'x'))
        self.assertTrue(hasattr(data, 'edge_index'))
        self.assertTrue(hasattr(data, 'pos'))

        self.assertEqual(data.x.ndim, 2)
        self.assertEqual(data.pos.ndim, 2)
        self.assertEqual(data.edge_index.shape[0], 2)
        self.assertEqual(data.x.shape[0], data.pos.shape[0])
        self.assertGreater(data.x.shape[0], 0)
        self.assertGreater(data.edge_index.shape[1], 0)

class TestResidueGraphBuilder(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.test_dir = os.path.dirname(__file__)
        cls.pdb_path = os.path.join(cls.test_dir, "pdbs", "3ux9.pdb")
        cls.cif_path = os.path.join(cls.test_dir, "pdbs", "3ux9.cif")
        cls.builder = ResidueGraphBuilder()

    def test_build_graph_from_pdb(self):
        data = self.builder.build(filename=self.pdb_path)
        self._validate_graph(data)

    def test_build_graph_from_cif(self):
        data = self.builder.build(filename=self.cif_path)
        self._validate_graph(data)

    def _validate_graph(self, data):
        self.assertIsInstance(data, Data)
        self.assertTrue(hasattr(data, 'x'))
        self.assertTrue(hasattr(data, 'edge_index'))
        self.assertTrue(hasattr(data, 'pos'))

        self.assertEqual(data.x.ndim, 2)
        self.assertEqual(data.pos.ndim, 2)
        self.assertEqual(data.edge_index.shape[0], 2)
        self.assertEqual(data.x.shape[0], data.pos.shape[0])
        self.assertGreater(data.x.shape[0], 0)
        self.assertGreater(data.edge_index.shape[1], 0)

if __name__ == "__main__":
    unittest.main()

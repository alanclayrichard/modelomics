import unittest
import torch
import os
from modelomics import sequence_embeddings as sequences

# Set TEST_MODEL env var to test specific models: "esmc", "e1", "esm3", "prott5", or "all"
# Example: TEST_MODEL=esmc python test_embeddings.py
TEST_MODEL = os.environ.get("TEST_MODEL", "esmc")

class TestEmbeddings(unittest.TestCase):
    def setUp(self):
        self.test_dir = os.path.dirname(__file__)
        self.cif_path = os.path.join(self.test_dir, "pdbs", "3ux9.cif")
        self.test_seq = "AAAAAAA"

    def test_sequence_from_cif(self):
        seq = sequences.sequence_from_cif(self.cif_path, chain='B')
        self.assertIsInstance(seq, str)
        self.assertGreater(len(seq), 0)

    def test_embed_sequence_esmc(self):
        embeddings = sequences.embed_sequence(self.test_seq)
        self.assertIsInstance(embeddings, torch.Tensor)
        self.assertEqual(len(embeddings.shape), 3)
        self.assertEqual(embeddings.shape[0], 1)
        self.assertEqual(embeddings.shape[1], len(self.test_seq)+2) 
        self.assertGreater(embeddings.shape[2], 0)

    def test_embed_sequence_esmc_masked(self):
        mask_sites = [1, 2]
        embeddings = sequences.embed_sequence(self.test_seq, mask_positions=mask_sites)
        self.assertIsInstance(embeddings, torch.Tensor)
        self.assertEqual(len(embeddings.shape), 3)
        self.assertEqual(embeddings.shape[0], 1)
        self.assertEqual(embeddings.shape[1], len(self.test_seq)+2)
        self.assertGreater(embeddings.shape[2], 0)

    def test_embed_sequence_e1(self):
        embeddings = sequences.embed_sequence(self.test_seq, model="e1")
        self.assertIsInstance(embeddings, torch.Tensor)
        self.assertEqual(len(embeddings.shape), 2)
        self.assertEqual(embeddings.shape[0], len(self.test_seq))
        self.assertGreater(embeddings.shape[1], 0)

    def test_embed_sequence_e1_300m(self):
        embeddings = sequences.embed_sequence(self.test_seq, model="e1", model_size="300m")
        self.assertIsInstance(embeddings, torch.Tensor)
        self.assertEqual(len(embeddings.shape), 2)
        self.assertEqual(embeddings.shape[0], len(self.test_seq))
        self.assertGreater(embeddings.shape[1], 0)

    def test_embed_sequence_esm3(self):
        embeddings = sequences.embed_sequence(self.test_seq, model="esm3")
        self.assertIsInstance(embeddings, torch.Tensor)
        self.assertEqual(len(embeddings.shape), 3)
        self.assertEqual(embeddings.shape[0], 1)
        self.assertEqual(embeddings.shape[1], len(self.test_seq)+2)
        self.assertGreater(embeddings.shape[2], 0)

    def test_embed_sequence_esm3_masked(self):
        mask_sites = [1, 2]
        embeddings = sequences.embed_sequence(self.test_seq, model="esm3", mask_positions=mask_sites)
        self.assertIsInstance(embeddings, torch.Tensor)
        self.assertEqual(len(embeddings.shape), 3)
        self.assertEqual(embeddings.shape[0], 1)
        self.assertEqual(embeddings.shape[1], len(self.test_seq)+2)
        self.assertGreater(embeddings.shape[2], 0)

    def test_embed_sequence_prott5(self):
        embeddings = sequences.embed_sequence(self.test_seq, model="prott5")
        self.assertIsInstance(embeddings, torch.Tensor)
        self.assertEqual(len(embeddings.shape), 2)
        self.assertEqual(embeddings.shape[0], len(self.test_seq))
        self.assertEqual(embeddings.shape[1], 1024)  # ProtT5 has 1024 dims

    def test_embed_sequence_invalid_model(self):
        with self.assertRaises(ValueError):
            sequences.embed_sequence(self.test_seq, model="invalid_model")

if __name__ == "__main__":
    unittest.main()
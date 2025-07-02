import unittest
import torch
import os
from modelomics import sequences

class TestEmbeddings(unittest.TestCase):
    def setUp(self):
        self.test_dir = os.path.dirname(__file__)
        self.cif_path = os.path.join(self.test_dir, "pdbs", "3ux9.cif")

    def test_sequence_from_cif(self):
        seq = sequences.sequence_from_cif(self.cif_path, chain='B')
        self.assertIsInstance(seq, str)
        self.assertGreater(len(seq), 0)

    def test_embed_sequences(self):
        seq = sequences.sequence_from_cif(self.cif_path, chain='B')
        embeddings = sequences.embed_sequence(seq)
        self.assertIsInstance(embeddings, torch.Tensor)
        self.assertEqual(len(embeddings.shape), 3)
        self.assertEqual(embeddings.shape[0], 1)
        self.assertEqual(embeddings.shape[1], len(seq)+2) 
        self.assertGreater(embeddings.shape[2], 0) 

    def test_embed_sequences_masked(self):
        seq = sequences.sequence_from_cif(self.cif_path, chain='B')
        mask_sites = [1, 2]
        embeddings = sequences.embed_sequence_with_mask(seq, mask_sites)
        self.assertIsInstance(embeddings, torch.Tensor)
        self.assertEqual(len(embeddings.shape), 3)
        self.assertEqual(embeddings.shape[0], 1)
        self.assertEqual(embeddings.shape[1], len(seq)+2)
        self.assertGreater(embeddings.shape[2], 0)

if __name__ == "__main__":
    unittest.main()
import unittest
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
        print(embeddings)

if __name__ == "__main__":
    unittest.main()
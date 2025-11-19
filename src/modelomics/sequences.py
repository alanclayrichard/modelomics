# created by clay 07/01/25

# Load the pretrained ESMC model when called (lazy imports for speed)
_client = None 

def get_client():
    from esm.models.esmc import ESMC
    global _client
    if _client is None:
        _client = ESMC.from_pretrained("esmc_300m").to("cpu")
    return _client

def sequence_from_cif(file_path, chain='A'):
    from Bio.PDB import MMCIFParser, Polypeptide
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", file_path)

    # Get the model (first model)
    model = structure[0]

    # Get the specified chain
    if chain not in model:
        raise ValueError(f"Chain {chain} not found in the structure.")
    chain = model[chain]

    # Extract polypeptide
    ppb = Polypeptide.PPBuilder()
    sequence = ""
    for pp in ppb.build_peptides(chain):
        sequence += pp.get_sequence()
    
    return str(sequence)

def sequence_from_pdb(file_path, chain='A'):
    from Bio.PDB import PDBParser, Polypeptide
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", file_path)

    # Get the model (first model)
    model = structure[0]

    # Get the specified chain
    if chain not in model:
        raise ValueError(f"Chain {chain} not found in the structure.")
    chain = model[chain]

    # Extract polypeptide
    ppb = Polypeptide.PPBuilder()
    sequence = ""
    for pp in ppb.build_peptides(chain):
        sequence += pp.get_sequence()
    
    return str(sequence)

def embed_sequence(seq):
    from esm.sdk.api import ESMProtein, LogitsConfig
    client = get_client()
    # Define the protein sequence
    protein = ESMProtein(sequence=(seq))

    # Encode the protein into ESMProteinTensor
    protein_tensor = client.encode(protein)

    # Get logits and embeddings
    logits_output = client.logits(
        protein_tensor, LogitsConfig(sequence=True, return_embeddings=True)
    )

    # Extract the embeddings tensor
    embeddings = logits_output.embeddings

    return embeddings

def embed_sequence_with_mask(seq, mask_positions):
    from esm.sdk.api import ESMProtein, LogitsConfig
    client = get_client()
    # define the sequence
    protein = ESMProtein(sequence=seq)
    # encode the protein sequence
    protein_tensor = client.encode(protein)

    # clone to avoid in-place edits
    sequence_tokens = protein_tensor.sequence.clone()  
    mask_token_id = client.tokenizer.mask_token_id

    for pos in mask_positions:
        if pos < 0 or pos >= len(sequence_tokens):
            raise IndexError(f"mask position {pos} is out of sequence bounds")
        # replace the sequence token with the mask token at pos
        sequence_tokens[pos] = mask_token_id

    # replace the encoded sequence with the masked sequence
    protein_tensor.sequence = sequence_tokens

    # look up the embeddings from the model
    logits_output = client.logits(
        protein_tensor, LogitsConfig(sequence=True, return_embeddings=True)
    )
    embeddings = logits_output.embeddings

    # return embeddings tensor
    return embeddings
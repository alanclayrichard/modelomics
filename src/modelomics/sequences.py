from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein, LogitsConfig
from Bio.PDB import MMCIFParser, Polypeptide

# Load the pretrained ESMC model
client = ESMC.from_pretrained("esmc_300m").to("cpu")

def sequence_from_cif(file_path, chain='A'):
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

def embed_sequence(seq):
    # Define the protein sequence
    protein = ESMProtein(
        sequence=(seq)
    )

    # Encode the protein into ESMProteinTensor
    protein_tensor = client.encode(protein)

    # Mask the first residue
    sequence_tokens = protein_tensor.sequence
    mask_token_id = client.tokenizer.mask_token_id  
    sequence_tokens[0] = mask_token_id 

    # Update the protein tensor with the masked sequence
    protein_tensor.sequence = sequence_tokens

    # Get logits and embeddings
    logits_output = client.logits(
        protein_tensor, LogitsConfig(sequence=True, return_embeddings=True)
    )

    # Extract the embeddings tensor
    embeddings = logits_output.embeddings

    return embeddings
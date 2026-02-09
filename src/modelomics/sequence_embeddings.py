# created by clay 07/01/25
'''
various sequence functionalities for use in protein tasks
'''

# lazy imports for speed
_esmc_client = None 
_esm3_client = None
_e1_models = {}
_prott5_model = None
_prott5_tokenizer = None


def _get_esmc_client():
    '''get the ESMC pretrained weights'''
    from esm.models.esmc import ESMC
    global _esmc_client
    if _esmc_client is None:
        _esmc_client = ESMC.from_pretrained("esmc_300m").to("cpu")
    return _esmc_client

def _get_esm3_model(model_size="esm3-sm-open-v1"):
    from esm.models.esm3 import ESM3
    global _esm3_client
    if _esm3_client is None:
        _esm3_client = ESM3.from_pretrained(model_size).to("cpu")
    return _esm3_client

def _get_e1_model(model_size="150m"):
    '''get the E1 pretrained weights'''
    from E1.modeling import E1ForMaskedLM
    global _e1_models
    if model_size not in _e1_models:
        _e1_models[model_size] = E1ForMaskedLM.from_pretrained(f"Profluent-Bio/E1-{model_size}").to("cpu")
        _e1_models[model_size].eval()
    return _e1_models[model_size]

def _get_prott5_model():
    '''get the ProtT5 pretrained weights'''
    try:
        from transformers import T5Tokenizer, T5EncoderModel
    except (ImportError, OSError) as e:
        raise ImportError(
            f"Failed to import transformers library for ProtT5: {e}. "
            "Install with: pip install transformers sentencepiece protobuf"
        )
    
    import torch
    global _prott5_model, _prott5_tokenizer
    if _prott5_model is None:
        device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
        _prott5_tokenizer = T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_half_uniref50-enc', do_lower_case=False)
        _prott5_model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc")
        _prott5_model = _prott5_model.to(device)
        # only GPUs support half-precision; use full-precision on CPU
        if device == torch.device("cpu"):
            _prott5_model = _prott5_model.to(torch.float32)
        _prott5_model.eval()
    return _prott5_model, _prott5_tokenizer

def sequence_from_cif(file_path, chain='A'):
    '''
    get the seuqnece from a cif file
    '''
    from Bio.PDB import MMCIFParser, Polypeptide
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", file_path)

    # get the model (first model)
    model = structure[0]

    # get the specified chain
    if chain not in model:
        raise ValueError(f"Chain {chain} not found in the structure.")
    chain = model[chain]

    # extract polypeptide
    ppb = Polypeptide.PPBuilder()
    sequence = ""
    for pp in ppb.build_peptides(chain):
        sequence += pp.get_sequence()
    
    return str(sequence)

def sequence_from_pdb(file_path, chain='A'):
    ''' 
    get the sequence from a pdb file
    '''
    from Bio.PDB import PDBParser, Polypeptide
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", file_path)

    # get the model (first model)
    model = structure[0]

    # get the specified chain
    if chain not in model:
        raise ValueError(f"Chain {chain} not found in the structure.")
    chain = model[chain]

    # extract polypeptide
    ppb = Polypeptide.PPBuilder()
    sequence = ""
    for pp in ppb.build_peptides(chain):
        sequence += pp.get_sequence()
    
    return str(sequence)

def embed_sequence(seq, model="esmc", model_size="150m", mask_positions=None):
    '''
    embed a string sequence using E1, ESMC, ESM3, or ProtT5 model
    
    args:
        seq: protein sequence string
        model: "e1", "esmc", "esm3", or "prott5" (default: "esmc")
        model_size: for E1 only - "150m" or "300m" (default: "150m")
        mask_positions: list of positions to mask (ESMC/ESM3 only)
    
    returns:
        embeddings tensor (L, E) for residues only
    '''
    import torch
    
    if model.lower() == "e1":
        from E1.batch_preparer import E1BatchPreparer
        e1_model = _get_e1_model(model_size)
        batch_preparer = E1BatchPreparer()
        batch = batch_preparer.get_batch_kwargs([seq], device="cpu")
        
        with torch.autocast("cpu", dtype=torch.bfloat16, enabled=False):
            outputs = e1_model(
                input_ids=batch["input_ids"],
                within_seq_position_ids=batch["within_seq_position_ids"],
                global_position_ids=batch["global_position_ids"],
                sequence_ids=batch["sequence_ids"],
                past_key_values=None,
                use_cache=False,
                output_attentions=False,
                output_hidden_states=False,
            )
        
        embeddings = outputs.embeddings
        residue_selector = ~(batch_preparer.get_boundary_token_mask(batch["input_ids"]))
        return embeddings[0, residue_selector[0]]
    
    elif model.lower() == "esmc":
        from esm.sdk.api import ESMProtein, LogitsConfig
        client = _get_esmc_client()
        protein = ESMProtein(sequence=seq)
        protein_tensor = client.encode(protein)
        
        if mask_positions is not None:
            sequence_tokens = protein_tensor.sequence.clone()
            mask_token_id = client.tokenizer.mask_token_id
            for pos in mask_positions:
                if pos < 0 or pos >= len(sequence_tokens):
                    raise IndexError(f"mask position {pos} is out of sequence bounds")
                sequence_tokens[pos] = mask_token_id
            protein_tensor.sequence = sequence_tokens
        
        logits_output = client.logits(
            protein_tensor, LogitsConfig(sequence=True, return_embeddings=True)
        )
        return logits_output.embeddings
    
    elif model.lower() == "esm3":
        from esm.sdk.api import ESMProtein, LogitsConfig
        client = _get_esm3_model()
        protein = ESMProtein(sequence=seq)
        protein_tensor = client.encode(protein)
        
        if mask_positions is not None:
            sequence_tokens = protein_tensor.sequence.clone()
            mask_token_id = client.tokenizer.mask_token_id
            for pos in mask_positions:
                if pos < 0 or pos >= len(sequence_tokens):
                    raise IndexError(f"mask position {pos} is out of sequence bounds")
                sequence_tokens[pos] = mask_token_id
            protein_tensor.sequence = sequence_tokens
        
        logits_output = client.logits(
            protein_tensor, LogitsConfig(sequence=True, return_embeddings=True)
        )
        return logits_output.embeddings
    
    elif model.lower() == "prott5":
        import re
        prott5_model, prott5_tokenizer = _get_prott5_model()
        device = next(prott5_model.parameters()).device
        
        # replace all rare/ambiguous amino acids by X
        seq_cleaned = seq.replace('U','X').replace('Z','X').replace('O','X').replace('B','X')
        seq_len = len(seq_cleaned)
        
        # introduce white-space between all amino acids
        processed_seq = " ".join(list(seq_cleaned))
        
        # tokenize sequence
        token_encoding = prott5_tokenizer([processed_seq], add_special_tokens=True, padding="longest")
        input_ids = torch.tensor(token_encoding['input_ids']).to(device)
        attention_mask = torch.tensor(token_encoding['attention_mask']).to(device)
        
        # generate embeddings
        try:
            with torch.no_grad():
                embedding_repr = prott5_model(input_ids=input_ids, attention_mask=attention_mask)
        except RuntimeError as e:
            raise RuntimeError(f"RuntimeError during embedding for sequence (L={seq_len}): {str(e)}")
        
        # extract residue embeddings: slice off special tokens, keep only residue embeddings
        # embedding shape is (batch_size, seq_len + special_tokens, embedding_dim)
        emb = embedding_repr.last_hidden_state[0, :seq_len]
        return emb.detach().cpu()
    
    else:
        raise ValueError(f"model must be 'e1', 'esmc', 'esm3', or 'prott5' got '{model}'")
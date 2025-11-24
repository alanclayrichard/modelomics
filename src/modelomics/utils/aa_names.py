# created by clay
''' 
amino acid names for 1->3 and 3->1 conversion
'''

aa_list = "ARNDCEQGHILKMFPSTWYV"

aa_3to1 = {
    "ALA": "A",
    "ARG": "R", 
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLU": "E",
    "GLN": "Q",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    # Non-standard amino acids
    "SEC": "U",  # Selenocysteine
    "PYL": "O",  # Pyrrolysine
    "ASX": "B",  # Asparagine or Aspartic acid
    "GLX": "Z",  # Glutamine or Glutamic acid
    "XAA": "X"   # Unknown or any amino acid
}

aa_1to3 = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN", 
    "D": "ASP",
    "C": "CYS",
    "E": "GLU",
    "Q": "GLN",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
    # Non-standard amino acids
    "U": "SEC",  # Selenocysteine
    "O": "PYL",  # Pyrrolysine
    "B": "ASX",  # Asparagine or Aspartic acid
    "Z": "GLX",  # Glutamine or Glutamic acid
    "X": "XAA"   # Unknown or any amino acid
}

# Standard 20 amino acids only
standard_aa_3to1 = {k: v for k, v in aa_3to1.items() 
                    if v in aa_list}
standard_aa_1to3 = {k: v for k, v in aa_1to3.items() 
                    if k in aa_list}

# Amino acid properties
charged_aa = {"R", "H", "K", "D", "E"}
positive_aa = {"R", "H", "K"}
negative_aa = {"D", "E"}
polar_aa = {"R", "N", "D", "C", "E", "Q", "H", "K", "S", "T", "Y"}
nonpolar_aa = {"A", "V", "I", "L", "M", "F", "W", "P", "G"}
aromatic_aa = {"F", "W", "Y", "H"}
aliphatic_aa = {"A", "V", "I", "L"}

def convert_3to1(three_letter):
    """
    convert three-letter amino acid code to one-letter code
    """
    return aa_3to1.get(three_letter.upper(), "X")

def convert_1to3(one_letter):
    """
    convert one-letter amino acid code to three-letter code
    """
    return aa_1to3.get(one_letter.upper(), "XAA")
from typing import Union


AMINOACID_ALPHABET_1TO3_MAP = {'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 
							'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile', 
							'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro', 
							'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'}
                            
MOLECULAR_MASS_MAP = {'A': 89.094, 'R': 174.203, 'N': 132.119, 'D': 133.104, 'C': 121.154, 
                  'E': 147.131, 'Q': 146.146, 'G': 75.067, 'H': 155.156, 'I': 131.175,
                  'L': 131.175, 'K': 146.189, 'M': 149.208, 'F': 165.192, 'P': 115.132,
                  'S': 105.093, 'T': 119.119, 'W': 204.228, 'Y': 181.191, 'V': 117.148}

AAS = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}

def is_prot(prot: str) -> bool:
    """
    Checks is given sequence a protein
    Arguments:
        prot (str) - aminoacid sequence of protein
    Return:
        bool if sequence is correct
        ValueError('Please  check proteins sequences') if there were wrong symbols
    """
    prot = prot.upper()
    uniq_AAS = set(prot)
    aa_test = (uniq_AAS <= set(AAS))
    if aa_test == 0:
        raise ValueError('Please  check proteins sequences')
    return True


def convert_1to3(prot: str) -> str:
    """
    Converts 1-symbol aminoacid sequence into 3-symbol aminoacid sequence.
    Arguments: 
        -prot (str) - aminoacid sequence in uppercase 1-symbol format
    Return: 
        -output (str) - aminoacid sequence in 3-symbol format.
    """
    output = ''
    if len(prot) > 0:
        for i in prot:
            if i in AMINOACID_ALPHABET_1TO3_MAP:
                output += AMINOACID_ALPHABET_1TO3_MAP[i]
            else:
                raise ValueError('Input format: aminoacids in uppercase 1-letter symbols')
    return output


def calculate_mm(prot: str) -> float:
    """
    Calculates molecular mass of protein.
    Argumets:
        -prot (str) - aminoacid sequence in uppercase 1-symbol format.
    Return:
        -output (float) - molecular mass in float format with 2 digits after dot.
    """
    prot_seq = set(prot)
    output = 0
    if len(prot) == 1:
        output = MOLECULAR_MASS_MAP[prot]
    else:
        for i in prot_seq:
            output += prot.count(i) * MOLECULAR_MASS_MAP[i]
    output -= 18.0153 * (len(prot) - 1)
    return round(output,3)


def count_aa_length(prot: str) -> int:
    """ 
    Counts the length of the sequence
     Arguments: 
      -prot (str) - the sequence, which length should be counted
     Return:  
      -int - the result of the count
    """
    return len(prot)


def count_nucl_length (prot: str) -> int: 
    """
    Counts the length of the nucleotide sequence that codes the inputted aminoacid sequence
     Arguments: 
      -prot (str) - the sequence, which coding nucleotide sequence length should be counted
     Return:
      -int - the result of the count
    """
    return len(prot) * 3


def count_aa_content(prot: str) -> dict:
    """
    Counts each aminoacid in protein and returns thire quantity

    Arguments: prot (str) - one of the input protein sequences was given by protein_tools
    Return: aa_content (dict) - dict of aminoacids and their quantity in protein
    """

    prot = prot.upper()
    AA_CONTENT = {'I': 0, 'T': 0, 'H': 0, 'S': 0, 'N': 0, 'K': 0, 'E': 0, 'D': 0, 'P': 0, 'F': 0,
 'R': 0, 'C': 0, 'V': 0, 'Q': 0, 'A': 0, 'Y': 0, 'W': 0, 'M': 0, 'G': 0, 'L': 0}
    for i in prot:
        if i in AAS:
            AA_CONTENT[i] += 1
    return AA_CONTENT


def count_extinction_280nm(prot: str) -> int:
    """
    Counts extinction in 280nm according to W, Y, C (cystine) number.

    Transforms prot sequence into dictionary using count_aa_content(prot) function.
    Uses the formula: e = 5500 * W + 1490 * Y + 125 * C
    Cystine number is counted roughly.

    Arguments: prot (str) - one of the input protein sequences
    Return: e (int) - result of counts: extinction coefficient at 280 nm

    """
    aa_cont_dict = count_aa_content(prot)

    w_number = aa_cont_dict.get('W')
    y_number = aa_cont_dict.get('Y')
    c_number = aa_cont_dict.get('C')

    e = 5500 * w_number + 1490 * y_number + 125*(c_number//2)
    return e


def protein_tools (function : str, *prots : str) -> Union[int, list, str]: 
    """
    Consists of several functions, is able to:
      -built in fuction 'is_prot': check whether the inputted sequence is a peptide 
      -'count_length': count the length of the sequence
      -'count_nucleotide_length': count the length of the coding nucleotide sequence of the inputted sequence
      -'count_MOLECULAR_MASS_MAP': count the molecular mass of the sequence
      -'convert_1_to_3': convert 1-letter input style into 3-letter and vice versa
      -'show_content': show the aminoacid content of the sequence
      -''count_extinction_280nm': Count extinction coefficient 280nm
     Arguments:
      -function (str) - the name of the action, the user wants to do on the sequence(s)
      -prots (str) - the sequence(s) that should be manipulated
     Return:
      -int - results of counts
      -list or str - result of convertation or showing the content

    """
    functions = {'count_length':count_aa_length, 'count_nucleotide_length':count_nucl_length,
                 'count_MOLECULAR_MASS_MAP':calculate_mm, 'show_content':count_aa_content, 'convert_1_to_3':convert_1to3,
                  'count_extinction_280nm':count_extinction_280nm }
    protein = []
    for prot in prots:
        is_prot(prot)
        protein.append(functions[function](prot))
    if len(protein) == 1:
        return protein[0]
    else:
        return protein
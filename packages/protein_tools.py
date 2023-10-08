from typing import Union
aminoacid_alphabet_1to3 = {'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 
							'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile', 
							'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro', 
							'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'}
                            
molecular_mass = {'A': 89.094, 'R': 174.203, 'N': 132.119, 'D': 133.104, 'C': 121.154, 
                  'E': 147.131, 'Q': 146.146, 'G': 75.067, 'H': 155.156, 'I': 131.175,
                  'L': 131.175, 'K': 146.189, 'M': 149.208, 'F': 165.192, 'P': 115.132,
                  'S': 105.093, 'T': 119.119, 'W': 204.228, 'Y': 181.191, 'V': 117.148}


def is_prot(prot: str) -> bool:
    """
    Checks is given sequence a protein
    Arguments:
        prot (str) - aminoacid sequence of protein
    Return:
        bool if sequence is correct
        ValueError('Please  check proteins sequences') if there were wrong symbols
    """
    aas = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}
    prot = prot.upper()
    uniq_aas = set(prot)
    aa_test = (uniq_aas <= aas)
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
            if i in aminoacid_alphabet_1to3:
                output += aminoacid_alphabet_1to3[i]
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
        output = molecular_mass[prot]
    else:
        for i in prot_seq:
            output += prot.count(i) * molecular_mass[i]
    output -= 18.0153*(len(prot)-1)
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
    return len(prot)*3


def count_aa_content(prot: str) -> dict:
    """
    Counts each aminoacid in protein and returns thire quantity

    Arguments: prot (str) - one of the input protein sequences was given by protein_tools
    Return: aa_content (dict) - dict of aminoacids and their quantity in protein
    """

    aas = 'ACDEFGHIKLMNPQRSTVWY'
    prot = prot.upper()
    aa_counter = [0] * 20
    for i in range(len(prot)):
        n = aas.index(prot[i])
        aa_counter[n] += 1

    aa_content = dict(zip(list(aas), aa_counter))
    return aa_content


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

    W_number = aa_cont_dict.get('W')
    Y_number = aa_cont_dict.get('Y')
    C_number = aa_cont_dict.get('C')

    if C_number == 0:
        e = 5500 * W_number + 1490 * Y_number
    else:
        e = 5500 * W_number + 1490 * Y_number + 125*(C_number//2)
    return e


def protein_tools (function : str, *prots : str) -> Union[int, list, str]: 
    """
    Consists of several functions, is able to:
      -built in fuction 'is_prot': check whether the inputted sequence is a peptide 
      -'count_length': count the length of the sequence
      -'count_nucleotide_length': count the length of the coding nucleotide sequence of the inputted sequence
      -'count_molecular_mass': count the molecular mass of the sequence
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
                 'count_molecular_mass':calculate_mm, 'show_content':count_aa_content, 'convert_1_to_3':convert_1to3,
                  'count_extinction_280nm':count_extinction_280nm }
    protein = []
    for prot in prots:
        is_prot(prot)
        protein.append(functions[function](prot))
    if len(protein) == 1:
        return protein[0]
    else:
        return protein
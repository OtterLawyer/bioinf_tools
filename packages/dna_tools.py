from typing import Union


DNA2RNA = {'A':'A', 'a':'a', 'T':'U', 't':'u', 'G':'G', 'g':'g', 'C':'C', 'c':'c'}
DNA2CDNA = {'A':'T', 'a':'t', 'T':'A', 't':'a', 'G':'C', 'g':'c', 'C':'G', 'c':'g'}
RNA2CRNA = {'A':'U', 'a':'u', 'U':'A', 'u':'a', 'G':'C', 'g':'c', 'C':'G', 'c':'g'}
COMMANDS = {'transcribe': transcribe, 'reverse': reverse, 'complement':complement, 'reverse_complement': reverse_complement, 'count_gc': count_gc, 
'dna_or_rna': dna_or_rna, 'get_sequence_length': get_sequence_length}


def transcribe(seq: str) -> str:
    """
    Transcribes dna seqeunce into rna sequence.
    Arguments: 
        -seq (str): dna seqeunce.
    Return: 
        -rna seqeunce in str format.
    """
    rna = ''
    for nucl in seq:
        rna += DNA2RNA[nucl]
    return rna
    
def reverse(seq: str) -> str:
    """
    Transforms dna/rna sequence into it's reverse.
    Arguments: 
        -seq (str): dna/rna seqeunce.
    Return: 
        -str - reverse dna/rna sequence.
    """
    return seq[::-1]

def complement(seq: str) ->str:
    """
    Transforms dna/rna sequence into it's complement.
    Arguments: 
        -seq (str): dna/rna seqeunce.
    Return: 
        -str - complement dna/rna sequence.
    """
    cstrand = ''
    if 'u' in seq.lower():
        nucl_complement_map = RNA2CRNA
    else:
        nucl_complement_map = DNA2CDNA
    for nucl in seq:
        cstrand += nucl_complement_map[nucl]
    return cstrand

def reverse_complement(seq: str) -> str:
    """
    Transforms dna/rna sequence into it's reverse complement.
    Arguments: 
        -seq (str): dna/rna seqeunce.
    Return: 
        -str - reverse complement dna/rna sequence.
    """
    reverse_cstrand = complement(reverse(seq))
    return reverse_cstrand
        

def count_gc(seq: str) -> float:
    """
    Counts gc content of sequence.
    Arguments: 
        -seq (str): dna/rna seqeunce.
    Return: 
        -float - percent of GC content in sequence.
    """
    gc = seq.lower().count('g') + seq.lower().count('c')
    answer = gc / len(seq) * 100
    return round(answer, 2)

def dna_or_rna(seq: str) -> str:
    """
    Determines whether the input is DNA or RNA.
    Arguments:
        -seq (str): dna/rna seqeunce.
    Return:
        -str - dna or rna based on content.
    """
    if set(seq.lower()) == set('atgc'):
        return 'dna'
    elif set(seq.lower()) == set('augc'):
        return 'rna'    

def get_sequence_length(seq: str) -> int:
    """
    Counts the length of the sequence
    Arguments:
        -seq (str): dna/rna seqeunce.
    Return:
        -int - length of the sequence.
    """
    return len(seq)

def check(args: list[str]):
    """
    Checks list of sequences for correct format. 
    Arguments:
        -args (list) - Input format: at least 1 sequence and command, command on last position of list,
        only DNA/RNA seqeunces.
    Return:
        -ValueError - if something does not fit arguments input format.
    """
    if len(args)  < 2:
        raise ValueError('Input sequences and command')
    elif args[-1] not in COMMANDS:
        raise ValueError('No such command')
    for seq in args[:-1]:
        if set(seq.lower()) == set('autgc') and set(seq.lower()) != set('atgc') and set(seq.lower()) != set('augc'):
            raise ValueError('Incorrect nucleic acid')
def run_dna_rna_tools(*args: list[str]) -> Union[str, float, list[str]]:
    """
    Consists of several functions, is able to:
      -check whether the inputted sequence is DNA or RNA
      -transcribe DNA sequence into RNA sequence
      -transform into reverse sequence
      -transform into complement sequence
      -transform into reverse complement sequence
      -count GC content of sequence
      -determine whether sequence DNA of RNA
      -count length of sequence
      -check if given sequence is RNA or DNA 
      -count the length of the sequence
     Arguments:
      -*args - list[str]: 
        -function (str) - at last position, the name of the action, the user wants to do on the sequence(s)
        -seqs (str) - up to last position, the sequence(s) that should be manipulated
     Return:
      -int - results of counts
      -list or str - result of convertation or showing the content
    
    """
    seqs = args[:-1]
    cmd = args[-1]
    output = []
    check(args)
    if len(seqs) == 1:
        return(COMMANDS[cmd](seqs[0]))
    else:
        for seq in seqs:
            output.append(COMMANDS[cmd](seq))
        return(output[0])
    
    
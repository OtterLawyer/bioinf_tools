from typing import Union
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
def count_qscore(score_string: str) -> float:
    """
    Counts quality score for inputed quality score string.
    Argumets:
        -score_string (str) - score string from fastq file.
    Return:
        -mean_quality (float) - mean quality of sequence based on quality score string.
    """
    total_score = 0
    for symb in set(score_string):
        total_score += score_string.count(symb) * (ord(symb) - 33)
    return total_score / len(score_string)

def fastq_tools(
        seqs: Dict[str, str], 
        gc_bounds: Union[int, tuple[int, int]] = (0, 100), 
        length_bounds: Union[int, tuple[int, int]] = (0, 2**32), 
        quality_treshholds: int = 0
) -> Dict[str, str]:
    """
    Filters fastq files by specifiable parametrs.
    Arguments:
        -seqs (dict[str,str]) - a dictionary consisting of fastq sequences. Key - string, sequence name. The value is a tuple of two strings: sequence and quality.
        -gc_bounds (int,tuple) - GC interval (in percent) for filtering, default is (0, 100). If input is single int - it will be ceiling (0, n).
        -length_bounds (int,tuple) - length interval for filtering. Works exactly as gc_bounds. Default is (0, 2**32)
        -quality_treshholds (int) - The threshold value of average read quality for the filter is 0 by default (phred33 scale).
    Return:
        -output (dict[str,str]) - filtered dictionary, which consists of entities that fulfill entered parametrs.
    """
    output = {}
    if type(gc_bounds) == int:
        lower_gc, upper_gc = 0, gc_bounds
    else:
        lower_gc, upper_gc = gc_bounds[0], gc_bounds[1]

    if type(length_bounds) == int:
        lower_len, upper_len = 0, length_bounds
    else:
        lower_len, upper_len = length_bounds[0], length_bounds[1]
    
    for name in seqs:
        seq, q = seqs[name]
        if not(lower_gc <= float(count_gc(seq)) <= upper_gc):
            break
        elif not(lower_len <= len(seq) <= upper_len):
            break
        elif not(quality_treshholds <= count_qscore(q)):
            break
        else:
            output[name] = seqs[name]
    return output
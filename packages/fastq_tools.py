from typing import Union
import os


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
        input_path: str,
        output_filename: str = None,
        gc_bounds: Union[int, tuple[int, int]] = (0, 100), 
        length_bounds: Union[int, tuple[int, int]] = (0, 2**32), 
        quality_treshholds: int = 0
        ) -> dict[str, str]:
    """
    Filters fastq files by specifiable parametrs.
    Arguments:
        -input_path(str) - path to input file.
        -output_filename (str) - name of output file, input file name will be used if none is given.
        -gc_bounds (int,tuple) - GC interval (in percent) for filtering, default is (0, 100). If input is single int - it will be ceiling (0, n).
        -length_bounds (int,tuple) - length interval for filtering. Works exactly as gc_bounds. Default is (0, 2**32)
        -quality_treshholds (int) - The threshold value of average read quality for the filter is 0 by default (phred33 scale).
    Return:
        -output (dict[str,str]) - filtered dictionary, which consists of entities that fulfill entered parametrs.
    """
    if output_filename is None:
        output_filename = input_path

    output = {}
    if type(gc_bounds) == int:
        lower_gc, upper_gc = 0, gc_bounds
    else:
        lower_gc, upper_gc = gc_bounds[0], gc_bounds[1]

    if type(length_bounds) == int:
        lower_len, upper_len = 0, length_bounds
    else:
        lower_len, upper_len = length_bounds[0], length_bounds[1]
    
    seqs = dict()
    with open(input_path) as f:
        line = f.readline().split('runid')[0]
        while line.startswith('@'):
            seq = f.readline()
            sep = f.readline()
            qual = f.readline()
            seqs[line] = (seq, qual)
            line = f.readline().split('runid')[0]
    for name in seqs:
        seq, q = seqs[name]
        if not(lower_gc <= float(count_gc(seq)) <= upper_gc):
            continue
        elif not(lower_len <= len(seq) <= upper_len):
            continue
        elif not(quality_treshholds <= count_qscore(q)):
            continue
        else:
            output[name] = seqs[name]
    return output, output_filename

def write_fastq(data):
    seqs = data[0]
    output_filename = data[1]
    OUTDIR = 'fastq_filtrator_results'
    
    if os.path.isdir(OUTDIR) is not True:
        os.mkdir(OUTDIR)

    if not output_filename.endswith('.fastq'):
        output_filename = out_filename + '.fastq'
    out_path = os.path.join(OUTDIR, output_filename)
    
    
    with open(out_path, 'w') as out_f:
        for name in seqs:
            out_f.write(name)
            out_f.write(seqs[name][0])
            out_f.write('+\n')
            out_f.write(seqs[name][1])
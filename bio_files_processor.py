def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = 'oneline_output.fasta') -> None:
    """
    Converts from multiline format to oneline in fasta file.
    Arguments:
    -input_fasta (str) - name of input fasta file.
    -output_fasta (str) - name of output fasta file. If none is given will write to oneline_output.fasta file
    Return: None
    """
    if output_fasta is None:
        output_fasta = input_fasta
    else:
        output_fasta += '.fasta'
    out = []
    names = []
    oneline = ''

    with open('sample.fasta') as f:
        for i in f:
            if i.startswith('>'):
                names.append(i)
                out.append(oneline)
                oneline = ''
                continue
            else:
                oneline += i.strip()
    out.append(oneline)
    
    with open(output_fasta, 'w') as f:
        for i in range(len(names)):
            f.write(names[i])
            out_fasta = out[i+1] + '\n'
            f.write(out_fasta)

def change_fasta_start_pos(input_fasta: str, shift: int = 1, output_fasta:str = None) -> None:
    """
    Takes in fasta file with one oneline genome and shifts starting position based on shift variable.
    Arguments:
    -input_fasta (str) - input file name.
    -shift (int) - new starting position.
    -output_fasta (str) - name of output file. If none is given will output file will have 'shift' before input_fasta name.
    Return: None
    """
    shift_genome = ''
    if output_fasta is None:
        output_fasta = 'shift_' + input_fasta
    else:
        output_fasta += '.fasta'
    with open(input_fasta) as f:
        name = f.readline()
        genome = f.readline().strip()
        shift_genome = genome[shift:] + genome[:shift]
    with open(output_fasta, 'w') as f:
        f.write(name)
        f.write(shift_genome)
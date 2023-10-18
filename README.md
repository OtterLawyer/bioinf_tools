bioinf_tools
==============
This is collection of tools for dna, rna, protein and fastq.

***table of contents:***
- [Introduction](#Introduction)
- [Installation](#Installation)
- [How to use](#Use)
- [Autor](#Autor)

# Introduction <a name="Introduction"></a>

**bioinf_tools** — is project with fuctions for dna, rna, protein and fastq. Documentation contais recomendations and examples for execution. You can simply copy and execute this project.
bioinf_tools consists of 3 tools:
- na_tools - tools for DNA and RNA.
- p_tools - tools for proteins.
- fq_tools - tools for filtering fastq files.
- w_fq - tool for writing output of fq_tools to fastq file.

# Installation <a name="Installation"></a>

To install this project you need to download bioinf_tools.py and packages folder.

# How to use <a name="Use"></a>

You need to import module bioinf_tools. In this example we will import bioinf_tools as bi for short.
```python
import bioinf_tools as bi
```
## na_tools
na_tools takes only DNA and RNA sequences and command to execute. It takes at least 2 arguments first: first to second to last arguments dna/rna sequences and last argument must be command. Arguments must be string type. Arguments other then DNA and RNA will be ignored.

If count of sequences is greater than 1 run_dna_rna_tools returns list, if count of seqeunces is 1 it returns str or int depending on command. 

***There are 7 functions:***
- transcribe

Returns string of DNA/RNA variant of inputed sequence.
- reverse

Returns reverse sequence of inputed DNA/RNA.
- complement

Returns compelement sequence of inputed DNA/RNA.
- reverse_complement

Returns reverse complement of inputed DNA/RNA.
- count_gc

Takes sequence as input in string format. Counts GC content and returns string of percent of GC content in sequence.
- dna_or_rna

Takes sequence as input in string format. Returns 'dna' or 'rna' string depending on content of input sequences. Precision is up to 2 digits after dot.
- get_sequence_length

Takes sequence as input in string format. Returns length of sequences as integer type.
```python
# We aldready imported bioinf_tools, if you did not uncomment next line
#import bioinf_tools as bi
bi.na_tools('ATG', 'transcribe') # 'AUG'
bi.na_tools('ATG', 'reverse') # 'GTA'
bi.na_tools('AtG', 'complement') # 'TaC'
bi.na_tools('ATg', 'reverse_complement') # 'cAT'
bi.na_tools('ATGC', 'count_gc') # '50%'
bi.na_tools('AtG', 'dna_or_rna') # 'dna'
bi.na_tools('ATG', 'aT', 'get_sequence_length') # ['3', '2']
```
## p_tools
Provide a tool with the sequence(s) of the protein(s) in 1-letter format (for example, DYKDDDDK) and the function needed. If you
occasionally write down a non-peptide sequence, the programm will return an error.  

Here is the catalogue of actions the user can choose: 

- *count_length*: gives the length(s) of the protein sequence(s)  
- *count_nucleotide_length*: counts the length(s) of the coding nucleotide sequence(s) of the protein sequence(s)  
- *count_molecular_mass*: calculates molecular mass of the input (the algorithm takes into consideration water mass and subtracts it)    
- *show_content*: shows the aminoacid content of the protein(s)  
- *convert_1_to_3*: converts 1-letter format into 3-letter one  
- *count_extinction_280nm*: counts the molar extinction coefficient (this function counts cystine contribution to extinction coefficient as two cysteins give 1 SS-bond) 
```python
bi.p_tools('count_length', 'DYKDDDDK') #8
bi.p_tools('count_nucleotide_length', 'DYKDDDDK') #24
bi.p_tools('count_molecular_mass', 'DYKDDDDK') #1012.982
bi.p_tools('show_content', 'DYKDDDDK') #{'A': 0, 'C': 0, 'D': 5, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 2, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 1}
bi.p_tools('convert_1_to_3', 'DYKDDDDK') #'AspTyrLysAspAspAspAspLys'
bi.p_tools('count_extinction_280nm', 'DYKDDDDK') #1490
```

## fastq_tools
Filters fastq files by specifiable parametrs.
### Arguments:
- input_path (str) - path to input fastq file.
- output_filename (str) - name of output file. Output file will be writen to fastq_filtrator_resuls folder, with name output_filename.fastq. If none name is given input name will be used.
- gc_bounds (int,tuple) - GC interval (in percent) for filtering, default is (0, 100). If input is single int - it will be ceiling (0, n).
- length_bounds (int,tuple) - length interval for filtering. Works exactly as gc_bounds. Default is (0, 2**32)
- quality_treshholds (int) - The threshold value of average read quality for the filter is 0 by default (phred33 scale).
### Return:
- output (dict[str,str]) - filtered dictionary, which consists of entities that fulfill entered parametrs.

## write_fastq
Writes output of fastq_tools to fastq file.
### Arguments:
- data (tuple(dict[str,str], str)) - output of fastq_tools. Dictionary where key - string, seqeunce name. 

As example we will use first read from 1_control_psbA3_2019_minq7.fastq (https://zenodo.org/record/3736457/files/1_control_psbA3_2019_minq7.fastq?download=1)
Which consist of 192 nucleotides, has 46.88% GC content and has mean quality score ~17.1.
```python
fastq = {'ee15a423-b008-44be-a4b2': ('GTTGTACTTCGTTCAATCGGTAGGTGTTTAACCGGATGGTCACGCCTACCGTGACAAAGAGATTGTCGGTGTCTTTGTGTTTCTGTTGGTGCTGATATTGCATTATGCATGAACGTAATGCCCATTAGTTGTGAATCCACCATGCGCGGAAGATAGAGCGACAGGCAAGTCACAAAGACACCGACAACTGTC', "##$&$&/035881()'$0&*('-.=;685()$.%($'%%&#&)+..0,&+&%.-/+,%&()$3:0&@09BF=>CC8(78029F7=<=)+@+.6CCFFC@-8%2579<B8;88412134,,;:8./,#1#&(%((09;B=??48<=<@79*-:B540,8=B=444:<571-B5=ED2.56;110.5+,*)%%*")}
bi.fq_tools(fastq_file, (50, 100),(100,200),16) #{}
bi.fq_tools(fastq_file, (40, 100),(195,200),16) #{}
bi.fq_tools(fastq_file, (40, 100),(100,200),18) #{}
bi.fq_tools(fastq_file, (40, 100),(100,200),16) #{'ee15a423-b008-44be-a4b2': ('GTTGTACTTCGTTCAATCGGTAGGTGTTTAACCGGATGGTCACGCCTACCGTGACAAAGAGATTGTCGGTGTCTTTGTGTTTCTGTTGGTGCTGATATTGCATTATGCATGAACGTAATGCCCATTAGTTGTGAATCCACCATGCGCGGAAGATAGAGCGACAGGCAAGTCACAAAGACACCGACAACTGTC',
  "##$&$&/035881()'$0&*('-.=;685()$.%($'%%&#&)+..0,&+&%.-/+,%&()$3:0&@09BF=>CC8(78029F7=<=)+@+.6CCFFC@-8%2579<B8;88412134,,;:8./,#1#&(%((09;B=??48<=<@79*-:B540,8=B=444:<571-B5=ED2.56;110.5+,*)%%*")}
bi.w_fq(bi.fq_tools('1_control_psbA3_2019_minq7.fastq', out))
# will create fastq_filtrator_resuls/out.fastq file with all reads.
```

# Autor <a name="Autor"></a>
Sivtsev Aleksei 
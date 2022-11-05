#!/usr/bin/env python
# Author: Rachel Ferina rferina@uoregon.edu

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
These functions will be useful on FASTA and FASTQ files.'''

__version__ = "0.8"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

from lib2to3.pytree import convert
from xmlrpc.client import boolean
import re

DNA_bases = "AGTCCCATGGCFCCNNA"
RNA_bases = "UUACGAGGUUAUAGNU"


def convert_phred(letter: str) -> int:
    """Inputs a single character, and returns it as a phred score."""
    return ord(letter) - 33


def qual_score(phred_score: str) -> float:
    """Takes in an unmodified string of phred quality
    scores. Returns the average quality score of the
    input phred string.
    """
    phred_sum = 0
    for letter in phred_score:
        # convert letter to phred score
        score = convert_phred(letter)
        phred_sum += score
    # calculate average
    avg_phred = phred_sum / len(phred_score)
    return avg_phred


def validate_base_seq(seq: str, RNAflag: bool=False) -> bool:
    '''Takes in a string of DNA or RNA. Returns True if string has
    only As, Gs, Cs, Ns, and Ts (or Us if RNAflag), otherwise returns False. Case insensitive.'''
    # make sequence uppercase
    seq = seq.upper()
    # if the sequence isn't RNA, count AGTC and see if equal to length of sequence
    if RNAflag == False:
        if seq.count('A') + seq.count('G') + seq.count('T') + seq.count('C') + seq.count('N') == len(seq):
            return True
    # if the sequence is RNA, count AGUC and see if equal to length of sequence
    else:
        if seq.count('A') + seq.count('G') + seq.count('U') + seq.count('C') + seq.count('N') == len(seq):
            return True
    return False


def gc_content(seq: str) -> float:
    '''Takes in a string (DNA). Returns GC content of the DNA sequence as a decimal between 0 and 1.'''
    # make sequence uppercase
    seq = seq.upper()
    # count Gs
    Gs = seq.count("G")
    # count Cs
    Cs = seq.count("C")
    return (Gs+Cs) / len(seq)


def oneline_fasta(file):
    '''Makes FASTA sequences on one line. Writes out to the file
    fa_one_line.fa. Returns the number of records so they can be
    manually compared to the number of header lines in the output file,
    to confirm the output file is accurate.'''
    # make dict with headers as keys and sequences as values
    seq_dict = {}
    with open(file, 'r') as fa:
        line_count = 0
        for line in fa:
            line_count +=1
            line = line.strip('\n')
            # only get header lines
            if line[0] == '>':
                header_line = line
            # populate dict with seq lines (non-header lines)
            else:    
                if header_line not in seq_dict:
                    seq_dict[header_line] = line
                else:
                    seq_dict[header_line] += line
    # write out to file
    fa_one_line = open('fa_one_line.fa', 'w')
    for keys,vals in seq_dict.items():
        fa_one_line.write(str(keys) + '\n' + str(vals) + '\n')
    fa_one_line.close()
    return len(seq_dict)


def reverse_complement(DNA_str: str) -> str:
    '''
    Takes in a string of a DNA sequence, and returns the reverse
    complement of the sequence in a new string. N's don't have
    a reverse complement, and remain the same.
    '''
    rev_str = ''
    comp_dict = {'G': 'C', 'C': 'G', 'A': 'T', 'T': 'A', 'N': 'N'}
    comp_list = []
    position = len(DNA_str)
    for nuc in range(len(DNA_str)):
        position -= 1
        rev_str += DNA_str[position]
    for base in rev_str:
        comp_list.append(comp_dict[base])
    rev_comp = ''.join(comp_list)
    return rev_comp


def strand_checker(bitwise):
    """
    Takes in SAM file bitwise flag. Returns strandedness (+ or -) if
    the read is mapped.
    """
    bitwise = int(bitwise)
    if ((bitwise & 16) == 16):
        rev_comp = True
        strandedness = 'negative'
    else:
        strandedness = 'positive' 
    return strandedness


def adjust_position(sam_position, strand, cigar):
    """
    Takes in the SAM file position, DNA strand, and CIGAR 
    string, parses the CIGAR string and returns the 5'
    start position.
    """
    # by default, five prime start is position from SAM file
    five_prime_start = int(sam_position)
    # if strand is positive, parse leftmost S for softclipping if there
    if (strand == 'positive'):
        left_s = re.findall(r"^(\d+)S", cigar)
        if left_s:
            five_prime_start -= int(left_s[0])
        else:
            return five_prime_start
    # if strand is negative, parse skipped region, deletion, match/mismatch, and rightmost S (softclipping)
    elif strand == 'negative':
        if 'N' in cigar:
            skip = re.findall(r"(\d+)N", cigar)
            skipped_list = [int(i) for i in skip]
            skipped = sum(skipped_list)
        else:
            skipped = 0
        if 'D' in cigar:
            deletion = re.findall(r"(\d+)D", cigar)
            deletion_list = [int(i) for i in deletion]
            deleted = sum(deletion_list)
        else:
            deleted = 0
        if 'M' in cigar:
            match = re.findall(r"(\d+)M", cigar)
            match_list = [int(i) for i in match]
            matched = sum(match_list)
        else:
            matched = 0
        if 'S' in cigar[-1]:
            right_s = re.findall(r"(\d+)S$", cigar)
            if right_s:
                right_s = int(right_s[0])
        else:
            right_s = 0
        five_prime_start += skipped + deleted + right_s + matched
    return five_prime_start

if __name__ == "__main__":
    # write tests for functions above
    assert convert_phred('A') == 32, "incorrect 'A' convert_phred score"
    assert convert_phred('F') == 37, "incorrect 'F' convert_phred score"
    assert convert_phred("@") == 31, "incorrect '@' convert_phred score"
    print('passed convert_phred tests')

    assert qual_score('HJIC2@JFFH$$') == 30.166666666666668, "qual_score produced incorrect average"
    print('passed qual_score test')

    assert validate_base_seq("ACTCGCCT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("ACNNGCNT") == True, "Validate base seq does not work on DNA with N"
    assert validate_base_seq("UACAUG", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("UANNNG", True) == True, "Validate base seq does not work on RNA with N"
    assert validate_base_seq("CTGUUA") == False, "Validate base seq worked on invalid sequence"
    assert validate_base_seq(DNA_bases) == False, 'Validate base seq does not work on DNA_bases'
    assert validate_base_seq(RNA_bases, True) == True, 'Validate base seq does not work on RNA_bases'
    print('validate_base_seq passed DNA and RNA tests')

    assert gc_content("GCGCCCG") == 1
    assert gc_content("TTTATAAA") == 0
    assert gc_content("GGCCTATA") == 0.5
    assert gc_content("GCTATAAT") == 0.25
    print("gc_content passed GC content tests")

    file_fa = oneline_fasta('/projects/bgmp/rferina/bioinfo/Bi621/PS/PS7/Danio_rerio.GRCz11.pep.all.fa')
    # assert path.exists('fa_one_line.fa') == True, "fa_one_line_.fa file not created"
    print('oneline_fasta passed tests')

    assert reverse_complement('GTAGCGTA') == 'TACGCTAC', 'Reverse complement is incorrect'
    assert reverse_complement('TGTTCCGT') == 'ACGGAACA', 'Reverse complement is incorrect'
    assert reverse_complement('AGAGTCCA') == 'TGGACTCT', 'Reverse complement is incorrect'
    print('reverse_complement passed tests')

    assert strand_checker(0) == 'positive', 'Strand checker is incorrect'
    assert strand_checker(107) == 'positive', 'Strand checker is incorrect'
    assert strand_checker(16) == 'negative', 'Strand checker is incorrect'
    assert strand_checker(83) == 'negative', 'Strand checker is incorrect'

    assert adjust_position(15, 'positive', '5S15M') == 10, "Adjust position is incorrect"
    assert adjust_position(100, 'positive', '20S59M') == 80, "Adjust position is incorrect"
    assert adjust_position(15, 'positive', '15M6D5S') == 15, "Adjust position is incorrect"
    assert adjust_position(10, 'negative', '3S17M') == 27, "Adjust position is incorrect"
    assert adjust_position(20, 'negative', '30M5D3S') == 58, "Adjust position is incorrect"
    assert adjust_position(10, 'negative', '4S10M5D50N5M5S') == 85, "Adjust position is incorrect"
    assert adjust_position(140, 'negative', '10M20N5D15N10S') == 200, "Adjust position is incorrect"
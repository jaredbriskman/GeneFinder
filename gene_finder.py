# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Jared Briskman

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq
startCodon = 'ATG'
stopCodon = ['TAA', 'TGA', 'TAG']
dna = load_seq("./data/X73525.fa")



def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('G')
    'C'
    >>> get_complement('T')
    'A'
    >>> get_complement('S')
    'Not a nucleotide'
    """

    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    else:
        return "Not a nucleotide"


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    complement = ''
    for char in dna:
        complement = get_complement(char) + complement #prepends the reverse
    return complement


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGAGA")
    'ATGAGA'
    """
    for i in xrange(len(dna)/3):
        if dna[i*3:i*3 + 3] in stopCodon:
            return dna[:i*3]
    return dna






def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGGAA")
    [ATGGAA]
    """
    n = 0
    listORFS = []
    while n < len(dna):
        if dna[n:n+3] == startCodon:
            listORFS.append(rest_of_ORF(dna[n:]))
            n += len(rest_of_ORF(dna[n:])) + 3
        else:
            n += 3
    return listORFS




def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    OrfList = []
    OrfList.extend(find_all_ORFs_oneframe(dna))
    OrfList.extend(find_all_ORFs_oneframe(dna[1:]))
    OrfList.extend(find_all_ORFs_oneframe(dna[2:]))
    return OrfList


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """

    bothStrandsList = []
    bothStrandsList.extend(find_all_ORFs(dna))
    bothStrandsList.extend(find_all_ORFs(get_reverse_complement(dna)))
    return bothStrandsList


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    return len(max(find_all_ORFs_both_strands(dna), key=len))


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF

        Can't doctest due to random nature.
        """
    maxLength = None
    for i in xrange(num_trials):
        localMax = longest_ORF(shuffle_string(dna))
        if localMax > maxLength:
            maxLength = localMax
    return maxLength

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    AA = ""
    for i in xrange(len(dna)/3):
        AA += aa_table[dna[i*3:i*3 + 3]]
    return AA

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.

        Can't doctest due to random nature.
    """
    threshold = longest_ORF_noncoding(dna, 1500)
    print threshold
    raw_list = find_all_ORFs_both_strands(dna)
    aa_list = [coding_strand_to_AA(x) for x in raw_list if len(x) > threshold]
    return aa_list

if __name__ == "__main__":
    import doctest
    doctest.run_docstring_examples(coding_strand_to_AA,globals())

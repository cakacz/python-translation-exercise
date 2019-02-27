#! /usr/bin/env python3

import sys

def translate_sequence(rna_sequence, genetic_code):
    """Translates a sequence of RNA into a sequence of amino acids.

    Translates `rna_sequence` into string of amino acids, according to the
    `genetic_code` given as a dict. Translation begins at the first position of
    the `rna_sequence` and continues until the first stop codon is encountered
    or the end of `rna_sequence` is reached.

    If `rna_sequence` is less than 3 bases long, or starts with a stop codon,
    an empty string is returned.
    """
    amino = ''
    if len(rna_sequence) < 3:
        return(amino)
    else:
        sequence=rna_sequence.upper()
        #As long as the sequence is at least three characters, will put the 
        #first three in one string and save sequence as the remaining characters
        #then looks up value in genetic_code
        while len(sequence) >= 3:
            codon = sequence[0:3]
            sequence = sequence[3:]
            #if the codon is a stop codon, returns current sequence and exits
            if genetic_code[codon] =='*':
                return(amino)
            else:
                 amino = amino + genetic_code[codon]
        return(amino)

def get_all_translations(rna_sequence, genetic_code):
    """Get a list of all amino acid sequences encoded by an RNA sequence.

    All three reading frames of `rna_sequence` are scanned from 'left' to
    'right', and the generation of a sequence of amino acids is started
    whenever the start codon 'AUG' is found. The `rna_sequence` is assumed to
    be in the correct orientation (i.e., no reverse and/or complement of the
    sequence is explored).

    The function returns a list of all possible amino acid sequences that
    are encoded by `rna_sequence`.

    If no amino acids can be translated from `rna_sequence`, an empty list is
    returned.
    """
    #slice off the first base and the first two bases to get alternate 
    #reading frames
    #onetrimmed = []
    #twotrimmed = []
    #threetrimmed = []
    seqone = []
    rf1 = rna_sequence
    rf2 = rna_sequence[1:]
    rf3 = rna_sequence[2:]
    #run all three reading frames through function to get sequence after start codon
    #while len(rf1) >= 3:
    onetrimmed = find_start_codon(rf1)
        #onetrimmed = onetrimmed.extend(rf1)
    #while len(rf2) >= 3:
    twotrimmed = find_start_codon(rf2)
        #onetrimmed = onetrimmed.extend(rf2)
    #while len(rf3) >= 3: 
    threetrimmed = find_start_codon(rf3)
        #onetrimmed = onetrimmed.extend(rf3) 
   #check for empty lists, lists with values move through if statement
   #empty lists move to else clause 
    #for i in onetrimmed:
        #new = [translate_sequence(onetrimmed, genetic_code)]
#        print(seqone)
    if onetrimmed:
        seqone = [translate_sequence(onetrimmed, genetic_code)]
    else:
        seqone = []
#        print(seqone)
    if twotrimmed:
        seqtwo = [translate_sequence(twotrimmed, genetic_code)]
#        print(seqtwo)
    else:
        seqtwo = []
#        print(seqtwo)
    if threetrimmed:
        seqthree = [translate_sequence(threetrimmed, genetic_code)]
#        print(seqthree)
    else:
        seqthree = []
#        print(seqthree)
    seqone.extend(seqtwo)
    seqone.extend(seqthree)
#        seqone.extend(new)
    return(seqone)

def find_start_codon(sequence):
    """Moves through a given sequence by threes until it finds 'AUG'.
    It will then slice off the start codon and return the sequence to be translated."""
    cond = 1
    none = []
    while cond == 1:
        codon = sequence[0:3]
        sequence = sequence[3:]
        #when it finds a start codon it cuts the codon and returns the remaining sequence
        if codon =='AUG':
            sequence = codon + sequence
            cond = 0
            return(sequence)
        #if the codon is not a start and there are fewer than three bases exits with empty list
        elif len(sequence) < 3:
            cond = 0
            return(none)
        else:
            pass

def get_reverse(sequence):
    """Reverse orientation of `sequence`.

    Returns a string with `sequence` in the reverse order.

    If `sequence` is empty, an empty string is returned.
    """
    #force upppercase
    sequence=sequence.upper()
    #if an empty string, will evaluate to false and move to else statement
    if sequence:
        #reverse the string by slicing
        seq_reverse= sequence[::-1]
        return(seq_reverse)
    else:
        return('')

def get_complement(sequence):
    """Get the complement of `sequence`.

    Returns a string with the complementary sequence of `sequence`.

    If `sequence` is empty, an empty string is returned.
    """
    #force uppercase
    sequence=sequence.upper()
    if sequence:
        #use translate method to change each base to its complement
        compa=sequence.replace('A', 'X')
        compu=compa.replace('U', 'A')
        comp_both=compu.replace('X', 'U')
        compg=comp_both.replace('G', 'Z')
        compc=compg.replace('C', 'G')
        complement=compc.replace('Z', 'C')
        return(complement)
    else:
        return('')

def reverse_and_complement(sequence):
    """Get the reversed and complemented form of `sequence`.

    Returns a string that is the reversed and complemented sequence
    of `sequence`.

    If `sequence` is empty, an empty string is returned.
    """
    #just run it through get_complement and get_reverse
    comp=get_complement(sequence)
    revcomp=get_reverse(comp)
    return(revcomp)


def get_longest_peptide(rna_sequence, genetic_code):
    """Get the longest peptide encoded by an RNA sequence.

    Explore six reading frames of `rna_sequence` (three reading frames of the
    current orientation, and the reversed and complemented form) and return (as
    a string) the longest sequence of amino acids that it encodes, according to
    the `genetic_code`.

    If no amino acids can be translated from `rna_sequence` nor its reverse and
    complement, an empty list is returned.
    """
    none = ''
    main = get_all_translations(rna_sequence, genetic_code)
    revcomp = reverse_and_complement(rna_sequence)
    alt = get_all_translations(revcomp, genetic_code)
    main.extend(alt)
    #longest = max(main, key = len)
    if main:
        return max(main, key = len)
    else:
        return(none)
    #maybe add list together?
    #now just find a way to compare lengths of each part and pull out the longest string
    #if all empty, just return the empty list

if __name__ == '__main__':
    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}
    rna_seq = ("AUG"
            "UAC"
            "UGG"
            "CAC"
            "GCU"
            "ACU"
            "GCU"
            "CCA"
            "UAU"
            "ACU"
            "CAC"
            "CAG"
            "AAU"
            "AUC"
            "AGU"
            "ACA"
            "GCG")
    longest_peptide = get_longest_peptide(rna_sequence = rna_seq,
            genetic_code = genetic_code)
    assert isinstance(longest_peptide, str), "Oops: the longest peptide is {0}, not a string".format(longest_peptide)
    message = "The longest peptide encoded by\n\t'{0}'\nis\n\t'{1}'\n".format(
            rna_seq,
            longest_peptide)
    sys.stdout.write(message)
    if longest_peptide == "MYWHATAPYTHQNISTA":
        sys.stdout.write("Indeed.\n")

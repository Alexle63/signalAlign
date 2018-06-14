#!/usr/bin/env python3
"""Functions to create a positions file."""
########################################################################
# File: makePositionsFiles.py
#  executable: makePositionsFiles.py
#
# Author: Andrew Bailey
# History: 5/21/18 Created
########################################################################

import re
import os
import string
import array
import subprocess
from collections import Counter
from signalalign.utils.parsers import read_fasta
from py3helpers.utils import find_substring_indices, all_string_permutations
from signalalign.utils import kmer_iterator, reverse_complement


def find_gatc_motifs(sequence):
    """Generate index of 'A' within the 'GATC' motifs in a nucleotide sequence

    :param sequence: since GATC motif is in DNA, expecting a DNA nucleotide sequence
    :return: generator yielding index of 'A' within the 'GATC'
    """
    return find_substring_indices(sequence.upper(), "GATC", offset=1)


def find_different_char_index(start_string, edit_string):
    """Compares standard and modified string and identifies the index of the modified character
        ex: find_char_difference_index("CCAGG","CFAGG") = 1

    :param start_string: starting string
    :param edit_string: string with a single character edit from start_string
    :return: index of string difference
    """
    assert len(start_string) == len(edit_string), ""
    pos = [i for i in range(len(start_string)) if start_string[i] != edit_string[i]]
    assert len(pos) == 1, "Only one character difference allowed. " \
                          "start_string={}, edit_string={}".format(start_string, edit_string)
    return pos[0]


def find_modification_index_and_character(canonical_motif, replacement_motif):
    """Compares canonical and modified motif and identifies the
            index of the modified nucleotide, canonical character, and replacement character.

    note. everything is converted to uppercase

    ex. find_modification_index_and_character("ATGC", "ETGC") = 0, "A", "E"

    :param canonical_motif: canonical nucleotide bases
    :param replacement_motif: replacement motif
    :return: mod_index, old_char, new_char
    """
    canonical_motif = canonical_motif.upper()
    replacement_motif = replacement_motif.upper()
    assert canonical_motif != replacement_motif, "Canonical motif cannot be the same as replacement motif"
    assert set(canonical_motif) <= set("ATGC"), "Canonical motif must only have canonical nucleotides"
    pos = find_different_char_index(canonical_motif, replacement_motif)
    old_char = canonical_motif[pos]
    new_char = replacement_motif[pos]
    return pos, old_char, new_char


# def make_positions_file(reference, output_path, *args):
#     """Creates a tsv file with the following format ("contig", "position", "strand", "change_from", "change_to").
#     Given a reference sequence and sets of sequence motif changes we report the location of each change.
#
#     ex: x : args = [("CCAGG","CFAGG"), ("CCTGG","CFTGG")]"""
#     with open(output_path, "w") as outfile:
#         for header, comment, sequence in read_fasta(reference):
#             for pair in args:
#                 motif = pair[0]
#                 modified = pair[1]
#                 # get pos, old character and the replacement character
#                 pos, old_char, new_char = find_modification_index_and_character(motif, modified)
#                 # get get rev_complement of motif and modified
#                 motif_comp = get_motif_REVcomplement(motif)
#                 # changed from rev complement to expand the alphabet and not contain
#                 # replacements to a single character, it can be different across motifs
#                 modified_comp = motif_comp[:rev_comp_pos] + new_char + \
#                                 motif_comp[rev_comp_pos+1:]
#
#                 seq_str_fwd = reference.replace(motif, modified)
#                 seq_str_bwd = reference.replace(motif_comp, modified_comp)
#                 nuc_positions = nuc_position(seq_str_fwd, new_char)
#                 for pos in nuc_positions:
#                     outfile.write(seq_name + "\t" + np.str(pos) + "\t" + "+" + "\t"
#                                   + old_char +"\t" + new_char + "\n")
#                 nuc_positions = nuc_position(seq_str_bwd, new_char)
#                 for pos in nuc_positions:
#                     outfile.write(seq_name + "\t" + np.str(pos) + "\t" + "-" + "\t"
#                                   + old_char +"\t" + new_char + "\n")

#
# def make_gatc_position_file(fasta, outfile):
#     outfile = os.path.abspath(outfile)
#     fH = open(outfile, 'w')
#     fH.write("X\t")
#
#     seq = get_first_seq(fasta)
#     for i in find_gatc_motifs(seq):
#         assert seq[i] == "A"
#         fH.write("{}\t".format(i))
#     fH.write("\n")
#     fH.write("X\t")
#     for i in find_gatc_motifs(seq):
#         t_pos = i + 1
#         assert seq[t_pos] == "T"
#         fH.write("{}\t".format(t_pos))  # + 1 because this is for the reverse complement
#     fH.write("\n")
#     fH.close()
#     return outfile

def replace_motifs_sequence_positions(sequence, motifs, overlap=False):
    """Edit nucleotide sequence using find and replace motifs

    note: we convert sequence to uppercase

    :param sequence: nucleotide sequence
    :param motifs: list of motif's which need to be replaced: eg [[find, replace]], [["CCAGG", "CEAGG"]]
    :param overlap: boolean option to look for motif overlaps
    """
    new_sequence = list(sequence)
    already_repaced_indexes = set()
    # gather motifs
    for motif_pair in motifs:
        assert len(motif_pair) == 2 and type(
            motif_pair) is list, "Motifs must be structured as list of lists, even for one motif find and replace"
        # find edit character and offset
        offset, old_char, substitution_char = find_modification_index_and_character(motif_pair[0], motif_pair[1])
        for index in find_substring_indices(sequence, motif_pair[0], offset=offset, overlap=overlap):
            new_sequence[index] = substitution_char
            # make sure that there is no overlapping assignments of characters
            assert index not in already_repaced_indexes, "Motifs has two different edits to a single nucleotide " \
                                                         "location. Check motifs {}".format(motifs)
            already_repaced_indexes.add(index)

    subst_sequence = ''.join(new_sequence).upper()
    return subst_sequence


def replace_periodic_sequence_positions(sequence, step_size, offset, substitution_char):
    """Edit nucleotide sequence using by replacing every 'step_size' nucleotides with an offset

    note: we convert sequence to uppercase

    eg: replace_periodic_sequence_positions("ATGCATGC", 3, 1, "F") = "AFGCFTGF"

    :param sequence: nucleotide sequence
    :param step_size: every 'step_size' locations the offset position is changed
    :param offset: the
    :param substitution_char: replacement character
    """
    assert offset < step_size, "Offset has to be less than step size"
    sequence = list(sequence)
    for i in range(offset, len(sequence), step_size):
        sequence[i] = substitution_char
    subst_sequence = ''.join(sequence).upper()

    return subst_sequence


def replace_periodic_reference_positions(reference_location, sub_fasta_path, step, offset, substitution_char='X'):
    """Edit and write a reference sequence to a specified path by replacing periodic characters

    note: if sub_fasta_path exists it will return the path without creating a new file

    :param reference_location: input reference
    :param sub_fasta_path: location of edited reference
    :param step: size of gap between substitution characters
    :param offset: offset of when to start creating substiutions
    :param substitution_char: character to replace original character
    """
    if os.path.isfile(sub_fasta_path):
        print("[substitute_reference_positions] Substituted reference fasta file exists: {}".format(
            sub_fasta_path))
        return sub_fasta_path
    else:
        print("[substitute_reference_positions] Creating substituted reference fasta file: {}".format(
            sub_fasta_path))
        # write
        with open(sub_fasta_path, 'w') as outfasta:
            for header, comment, sequence in read_fasta(reference_location):
                subst_sequence = replace_periodic_sequence_positions(sequence, step, offset, substitution_char)
                print(
                    ">%s %s\n%s" % (header, "substituted:{},step:{},offset:{}".format(substitution_char, step, offset),
                                    subst_sequence), file=outfasta)

    return sub_fasta_path


def replace_motif_reference_positions(reference_location, sub_fasta_path, motifs, overlap=False):
    """Replace motif  reference sequence to a specific path

    :param reference_location: input reference
    :param sub_fasta_path: location of edited reference
    :param motifs: list of motif's which need to be replaced: eg [[find, replace]], [["CCAGG", "CEAGG"]]
    :param overlap: of overlap is possible, replace with overlap: eg [["AAA", "AAT"]] :  AAAA -> AATT
    """
    if os.path.isfile(sub_fasta_path):
        print("[substitute_reference_positions] Substituted reference fasta file exists: {}".format(
            sub_fasta_path))
        return sub_fasta_path
    else:
        print("[substitute_reference_positions] Creating substituted reference fasta file: {}".format(
            sub_fasta_path))
        # write
        with open(sub_fasta_path, 'w') as outfasta:
            for header, comment, sequence in read_fasta(reference_location):
                subst_sequence = replace_motifs_sequence_positions(sequence, motifs, overlap)
                print(">%s %s\n%s" % (header, "substituted:{}".format(motifs),
                                      subst_sequence), file=outfasta)
    return sub_fasta_path


def samtools_faidx_fasta(fasta_path, log=None):
    """Index fasta using samtools faidx

    note: samtools must be in PATH

    :param fasta_path: path to fasta file
    """
    # produce faidx file
    assert os.path.isfile(fasta_path), "Path to fasta file does not exist"
    index_path = "{}.fai".format(fasta_path)
    if not os.path.exists(index_path):
        if log:
            print("[{}] indexing reference {}".format(log, fasta_path))
        args = ["samtools", "faidx", fasta_path]
        subprocess.check_call(args)
    assert os.path.isfile(index_path), "Error creating FAIDX file for: {}".format(fasta_path)
    return index_path


def count_all_sequence_kmers(seq, k=5, rev_comp=False):
    """Count all the 5'-3' kmers of a nucleotide sequence, rev_comp counts rev_comp seq IN ADDITION to given sequence

    :param seq: nucleotide sequence
    :param k: size of kmer
    :param rev_comp: boolean option to count reverse complement kmers as well
    :return: dictionary of kmers with counts as values
    """
    # loop through kmers
    kmers = Counter()
    for kmer in kmer_iterator(seq, k):
        kmers[kmer] += 1
    if rev_comp:
        # loop through rev_comp kmers
        seq1 = reverse_complement(seq, reverse=True, complement=True)
        for kmer in kmer_iterator(seq1, k):
            kmers[kmer] += 1

    return kmers


def get_sequence_kmers(seq, k=5, rev_comp=False):
    """Get the set of all kmers from a sequence.

    :param seq: nucleotide sequence
    :param k: size of kmer
    :param rev_comp: boolean option to count reverse complement kmers as well
    :return: set of kmers
    """
    return set(count_all_sequence_kmers(seq, k=k, rev_comp=rev_comp).keys())


def get_motif_kmers(motif_pair, k, alphabet="ATGC"):
    """Given a motif pair, create a list of all kmers which contain modification"""
    assert len(motif_pair) == 2, "Motif pair must be a list of length 2. len(motif_pair) = {}".format(len(motif_pair))
    canonical = motif_pair[0]
    modified = motif_pair[1]
    motif_len = len(canonical)
    # get mod index and chars
    mod_index, old_char, new_char = find_modification_index_and_character(canonical, modified)
    bases_after = motif_len - mod_index - 1

    # get overlaps for front and back of kmer
    front_overlap, back_overlap = get_front_back_kmer_overlap(k, motif_len, mod_index)
    # pre-compute kmers
    kmer_set_dict = dict()
    for i in range(1, max(front_overlap, back_overlap)+1):
        kmer_set_dict[i] = [x for x in all_string_permutations(alphabet, i)]
    kmer_set_dict[0] = ['']

    motif_kmers = []
    for i in range(k):
        # get prepend kmers and index for front of motif
        if i >= front_overlap:
            front_index = i - front_overlap
            prepend_kmers = ['']
        else:
            prepend_kmers = kmer_set_dict[front_overlap-i]
            front_index = 0
        # get append kmers and index for back of motif
        if i > bases_after:
            append_kmers = kmer_set_dict[i - bases_after]
            back_index = motif_len
        else:
            back_index = mod_index+i+1
            append_kmers = ['']

        kmer = modified[front_index:back_index]
        motif_kmers.extend([front+kmer+back for front in prepend_kmers for back in append_kmers if front+kmer+back is not ''])

    return set(motif_kmers)


def get_front_back_kmer_overlap(k, motif_len, mod_index):
    """Get the largest number of bases overlap at front and back of motif

    eg: k=3 , motif_len = 2, mod_index = 1
        motif = GE

        X G E X
        _ _ _       max front_overlap = 1
          _ _ _     max back_overlap = 1

    :param k: length of kmer
    :param motif_len: length of motif
    :param mod_index: index position of modification
    :return: largest overlap in the front and back of a generated kmer
    """
    assert k >= 1, "k cannot be less than 1. k: {}".format(k)
    front_overlap = k - mod_index - 1
    back_overlap = k - (motif_len - mod_index)
    return front_overlap, back_overlap


# def wrapper(func, *args, **kwargs):
#     def wrapped():
#         return func(*args, **kwargs)
#     return wrapped


# def get_methyl_char(degenerate):
#     if degenerate == "adenosine":
#         return "I"
#     else:
#         return "E"
#
#
# def make_positions_file(fasta, degenerate, outfile):
#     if degenerate == "adenosine":
#         return make_gatc_position_file(fasta, outfile)
#     else:
#         return make_CCWGG_positions_file(fasta, outfile)
#
#
# def gatc_kmers(sequence_kmers, kmerlength):
#     assert kmerlength == 5 or kmerlength == 6, "only works with kmer lengths 5 and 6"
#     # NNNNGATCNNN
#     methyl_core = "GITC"
#     normal_core = "GATC"
#     nucleotides = "ACGT"
#
#     fourmers = [''.join(x) for x in product(nucleotides, repeat=4)]
#     threemers = [''.join(x) for x in product(nucleotides, repeat=3)]
#     twomers = [''.join(x) for x in product(nucleotides, repeat=2)]
#
#     labeled_kmers = []
#
#     # add NNNNGA*
#     if kmerlength == 6:
#         for fourmer in fourmers:
#             labeled_kmer = (fourmer + methyl_core)[:kmerlength]
#             normal_kmer = (fourmer + normal_core)[:kmerlength]
#             if normal_kmer in sequence_kmers:
#                 labeled_kmers.append(labeled_kmer)
#     # add NNNGA*T and NNNGA*
#     for threemer in threemers:
#         labeled_kmer = (threemer + methyl_core)[:kmerlength]
#         normal_kmer = (threemer + normal_core)[:kmerlength]
#         if normal_kmer in sequence_kmers:
#             labeled_kmers.append(labeled_kmer)
#         # A*TCNNN
#         if kmerlength == 6:
#             labeled_kmer = (methyl_core + threemer)[1:]
#             normal_kmer = (normal_core + threemer)[1:]
#             if normal_kmer in sequence_kmers:
#                 labeled_kmers.append(labeled_kmer)
#     # add NNGA*TC and NNGA*T
#     for twomer in twomers:
#         labeled_kmer = (twomer + methyl_core)[:kmerlength]
#         normal_kmer = (twomer + normal_core)[:kmerlength]
#         if normal_kmer in sequence_kmers:
#             labeled_kmers.append(labeled_kmer)
#         # A*TCNN
#         if kmerlength == 5:
#             labeled_kmer = (methyl_core + twomer)[1:]
#             normal_kmer = (normal_core + twomer)[1:]
#             if normal_kmer in sequence_kmers:
#                 labeled_kmers.append(labeled_kmer)
#         # NGA*TCN
#         if kmerlength == 6:
#             labeled_kmer = (twomer[0] + methyl_core + twomer[1])
#             normal_kmer = (twomer[0] + normal_core + twomer[1])
#             if normal_kmer in sequence_kmers:
#                 labeled_kmers.append(labeled_kmer)
#             # TODO forgot GA*TCNN for 6mers!
#             labeled_kmer = (methyl_core + twomer)
#             normal_kmer = (normal_core + twomer)
#             if normal_kmer in sequence_kmers:
#                 labeled_kmers.append(labeled_kmer)
#     if kmerlength == 5:
#         for onemer in "ACTG":
#             labeled_kmer = onemer + methyl_core
#             normal_kmer = onemer + normal_core
#             if normal_kmer in sequence_kmers:
#                 labeled_kmers.append(labeled_kmer)
#             labeled_kmer = methyl_core + onemer
#             normal_kmer = normal_core + onemer
#             if normal_kmer in sequence_kmers:
#                 labeled_kmers.append(labeled_kmer)
#     return set(labeled_kmers)
#
#
# def ctag_kmers(sequence_kmers, kmerlength):
#     assert kmerlength == 5 or kmerlength == 6, "only works with kmer lengths 5 and 6"
#     # NNNCTAGNNNN
#     methyl_core = "CTIG"
#     normal_core = "CTAG"
#     nucleotides = "ACGT"
#
#     fourmers = [''.join(x) for x in product(nucleotides, repeat=4)]
#     threemers = [''.join(x) for x in product(nucleotides, repeat=3)]
#     twomers = [''.join(x) for x in product(nucleotides, repeat=2)]
#
#     labeled_kmers = []
#
#     # add A*GNNNN
#     if kmerlength == 6:
#         for fourmer in fourmers:
#             labeled_kmer = (methyl_core + fourmer)[2:]
#             normal_kmer = (normal_core + fourmer)[2:]
#             if normal_kmer in sequence_kmers:
#                 labeled_kmers.append(labeled_kmer)
#     # add NNNCTA*
#     for threemer in threemers:
#         if kmerlength == 6:
#             labeled_kmer = (threemer + methyl_core)[:kmerlength]
#             normal_kmer = (threemer + normal_core)[:kmerlength]
#             if normal_kmer in sequence_kmers:
#                 labeled_kmers.append(labeled_kmer)
#             labeled_kmer = (methyl_core + threemer)[1:]
#             normal_kmer = (normal_core + threemer)[1:]
#             if normal_kmer in sequence_kmers:
#                 labeled_kmers.append(labeled_kmer)
#         # A*GNNN
#         if kmerlength == 5:
#             labeled_kmer = (methyl_core + threemer)[2:]
#             normal_kmer = (normal_core + threemer)[2:]
#             if normal_kmer in sequence_kmers:
#                 labeled_kmers.append(labeled_kmer)
#
#     # add NNCTA*G and NNCTA*
#     for twomer in twomers:
#         labeled_kmer = (twomer + methyl_core)[:kmerlength]
#         normal_kmer = (twomer + normal_core)[:kmerlength]
#         if normal_kmer in sequence_kmers:
#             labeled_kmers.append(labeled_kmer)
#         # CTA*GNN
#         if kmerlength == 6:
#             labeled_kmer = (methyl_core + twomer)[:kmerlength]
#             normal_kmer = (normal_core + twomer)[:kmerlength]
#             if normal_kmer in sequence_kmers:
#                 labeled_kmers.append(labeled_kmer)
#         # TA*GNN
#         if kmerlength == 5:
#             labeled_kmer = (methyl_core + twomer)[1:]
#             normal_kmer = (normal_core + twomer)[1:]
#             if normal_kmer in sequence_kmers:
#                 labeled_kmers.append(labeled_kmer)
#
#     if kmerlength == 5:
#         for onemer in nucleotides:
#             labeled_kmer = onemer + methyl_core
#             normal_kmer = onemer + normal_core
#             if normal_kmer in sequence_kmers:
#                 labeled_kmers.append(labeled_kmer)
#             labeled_kmer = methyl_core + onemer
#             normal_kmer = normal_core + onemer
#             if normal_kmer in sequence_kmers:
#                 labeled_kmers.append(labeled_kmer)
#     return set(labeled_kmers)
#
#
# def ccwgg_kmers(sequence_kmers, kmer_length):
#     def check_and_add(methyl_kmer):
#         normal_kmer = string.translate(methyl_kmer, demethylate)
#         if normal_kmer in sequence_kmers:
#             labeled_kmers.append(methyl_kmer)
#
#     labeled_kmers = []
#
#     methyl_core1 = "CEAGG"
#     methyl_core2 = "CETGG"
#     demethylate = string.maketrans("E", "C")
#
#     nucleotides = "ACGT"
#     fourmers = [''.join(x) for x in product(nucleotides, repeat=4)]
#     threemers = [''.join(x) for x in product(nucleotides, repeat=3)]
#     twomers = [''.join(x) for x in product(nucleotides, repeat=2)]
#     # NNNNCC*WGGNN
#
#     # NNNNCC*
#     if kmer_length == 6:
#         for fourmer in fourmers:
#             labeled_kmer1 = (fourmer + methyl_core1)[:kmer_length]
#             labeled_kmer2 = (fourmer + methyl_core2)[:kmer_length]
#             check_and_add(labeled_kmer1)
#             check_and_add(labeled_kmer2)
#
#     # NNNCC*W and NNNCC*
#     for threemer in threemers:
#         labeled_kmer1 = (threemer + methyl_core1)[:kmer_length]
#         labeled_kmer2 = (threemer + methyl_core2)[:kmer_length]
#         check_and_add(labeled_kmer1)
#         check_and_add(labeled_kmer2)
#
#     # NNCC*WG and NNCC*W
#     for twomer in twomers:
#         labeled_kmer1 = (twomer + methyl_core1)[:kmer_length]
#         labeled_kmer2 = (twomer + methyl_core2)[:kmer_length]
#         check_and_add(labeled_kmer1)
#         check_and_add(labeled_kmer2)
#         # C*WGGNN
#         if kmer_length == 6:
#             labeled_kmer1 = (methyl_core1 + twomer)[1:]
#             labeled_kmer2 = (methyl_core2 + twomer)[1:]
#             check_and_add(labeled_kmer1)
#             check_and_add(labeled_kmer2)
#
#     for onemer in nucleotides:
#         # CC*WGGN and C*WGGN
#         labeled_kmer1 = methyl_core1 + onemer
#         labeled_kmer2 = methyl_core2 + onemer
#         if kmer_length == 6:
#             check_and_add(labeled_kmer1)
#             check_and_add(labeled_kmer2)
#         if kmer_length == 5:
#             check_and_add(labeled_kmer1[1:])
#             check_and_add(labeled_kmer2[1:])
#         labeled_kmer1 = (onemer + methyl_core1)[:kmer_length]
#         labeled_kmer2 = (onemer + methyl_core2)[:kmer_length]
#         check_and_add(labeled_kmer1)
#         check_and_add(labeled_kmer2)
#
#     if kmer_length == 5:
#         check_and_add(methyl_core1)
#         check_and_add(methyl_core2)
#
#     return set(labeled_kmers)
#
#
# def ggwcc_kmers(sequence_kmers, kmer_length):
#     def check_and_add(methyl_kmer):
#         normal_kmer = string.translate(methyl_kmer, demethylate)
#         if normal_kmer in sequence_kmers:
#             labeled_kmers.append(methyl_kmer)
#
#     labeled_kmers = []
#
#     methyl_core1 = "GGAEC"
#     methyl_core2 = "GGTEC"
#     demethylate = string.maketrans("E", "C")
#
#     nucleotides = "ACGT"
#     fourmers = [''.join(x) for x in product(nucleotides, repeat=4)]
#     threemers = [''.join(x) for x in product(nucleotides, repeat=3)]
#     twomers = [''.join(x) for x in product(nucleotides, repeat=2)]
#
#     # NNGGWC*CNNN
#
#     # C*CNNNN
#     for fourmer in fourmers:
#         labeled_kmer1 = (methyl_core1 + fourmer)[3:]
#         labeled_kmer2 = (methyl_core2 + fourmer)[3:]
#         check_and_add(labeled_kmer1)
#         check_and_add(labeled_kmer2)
#
#     # WC*CNNN and C*CNNN
#     for threemer in threemers:
#         labeled_kmer1 = (methyl_core1 + threemer)[2:] if kmer_length == 6 else (methyl_core1 + threemer)[3:]
#         labeled_kmer2 = (methyl_core2 + threemer)[2:] if kmer_length == 6 else (methyl_core2 + threemer)[3:]
#         check_and_add(labeled_kmer1)
#         check_and_add(labeled_kmer2)
#
#     # GWC*CNN and WC*CNN
#     for twomer in twomers:
#         labeled_kmer1 = (methyl_core1 + twomer)[1:] if kmer_length == 6 else (methyl_core1 + twomer)[2:]
#         labeled_kmer2 = (methyl_core2 + twomer)[1:] if kmer_length == 6 else (methyl_core2 + twomer)[2:]
#         check_and_add(labeled_kmer1)
#         check_and_add(labeled_kmer2)
#         # NNGGWC*
#         if kmer_length == 6:
#             labeled_kmer1 = (twomer + methyl_core1)[:kmer_length]
#             labeled_kmer2 = (twomer + methyl_core2)[:kmer_length]
#             check_and_add(labeled_kmer1)
#             check_and_add(labeled_kmer2)
#
#     for onemer in nucleotides:
#         # NGGWC* and NGGWC*C
#         labeled_kmer1 = (onemer + methyl_core1)[:kmer_length]
#         labeled_kmer2 = (onemer + methyl_core2)[:kmer_length]
#         check_and_add(labeled_kmer1)
#         check_and_add(labeled_kmer2)
#         # GGWC*CN GWC*CN
#         labeled_kmer1 = methyl_core1 + onemer if kmer_length == 6 else (methyl_core1 + onemer)[1:]
#         labeled_kmer2 = methyl_core2 + onemer if kmer_length == 6 else (methyl_core2 + onemer)[1:]
#         check_and_add(labeled_kmer1)
#         check_and_add(labeled_kmer2)
#
#     if kmer_length == 5:
#         check_and_add(methyl_core1)
#         check_and_add(methyl_core2)
#
#     return set(labeled_kmers)
#
#
# def motif_kmers(core, kmer_length=5, multiplier=5):
#     motifs = []
#     repeat = kmer_length - len(core)
#     if repeat == 0:
#         return [core] * multiplier
#     else:
#         for k in product("ACGT", repeat=repeat):
#             fix = ''.join(k)
#             motifs.append(fix + core)
#             motifs.append(core + fix)
#         return motifs * multiplier
#
#


#
#
# def make_gatc_or_ccwgg_motif_file(fasta, degenerate, outfile):
#     if degenerate == "adenosine":
#         motif_finder = find_gatc_motifs
#     else:
#         motif_finder = find_ccwgg_motifs
#     outfile = os.path.abspath(outfile)
#     seq = get_first_seq(fasta)
#     positions = [x for x in motif_finder(seq)]
#     make_motif_file(positions, seq, outfile)
#     return outfile
#

# TODO write these tests ya dig
def getFastaDictionary(fastaFile):
    """Returns a dictionary of the first words of fasta headers to their corresponding
    fasta sequence
    """
    namesAndSequences = [(x[0].split()[0], x[1]) for x in fastaRead(open(fastaFile, 'r'))]
    names = [x[0] for x in namesAndSequences]
    assert len(names) == len(set(names))  # Check all the names are unique
    return dict(namesAndSequences)  # Hash of names to sequences


def fastaRead(fileHandleOrFile):
    """iteratively yields a sequence for each '>' it encounters, ignores '#' lines
    """
    fileHandle = _getFileHandle(fileHandleOrFile)
    line = fileHandle.readline()
    chars_to_remove = "\n "
    valid_chars = {x for x in string.ascii_letters + "-"}
    while line != '':
        if line[0] == '>':
            name = line[1:-1]
            line = fileHandle.readline()
            seq = array.array('b')
            while line != '' and line[0] != '>':
                line = line.translate(str.maketrans('', '', chars_to_remove))
                if len(line) > 0 and line[0] != '#':
                    seq.extend(list(map(ord, line)))
                line = fileHandle.readline()
            try:
                assert all(chr(x) in valid_chars for x in seq)
            except AssertionError:
                bad_chars = {chr(x) for x in seq if chr(x) not in valid_chars}
                raise RuntimeError("Invalid FASTA character(s) see in fasta sequence: {}".format(bad_chars))
            yield name, seq.tobytes()
        else:
            line = fileHandle.readline()
    if isinstance(fileHandleOrFile, "".__class__):
        fileHandle.close()


def _getFileHandle(fileHandleOrFile, mode="r"):
    if isinstance(fileHandleOrFile, "".__class__):
        return open(fileHandleOrFile, mode)
    else:
        return fileHandleOrFile


def fastaWrite(fileHandleOrFile, name, seq, mode="w"):
    """Writes out fasta file
    """
    fileHandle = _getFileHandle(fileHandleOrFile, mode)
    valid_chars = {x for x in string.ascii_letters + "-"}
    try:
        assert any([isinstance(seq, str), isinstance(seq, str)])
    except AssertionError:
        raise RuntimeError("Sequence is not unicode or string")
    try:
        assert all(x in valid_chars for x in seq)
    except AssertionError:
        bad_chars = {x for x in seq if x not in valid_chars}
        raise RuntimeError("Invalid FASTA character(s) see in fasta sequence: {}".format(bad_chars))
    fileHandle.write(">%s\n" % name)
    chunkSize = 100
    for i in range(0, len(seq), chunkSize):
        fileHandle.write("%s\n" % seq[i:i + chunkSize])
    if isinstance(fileHandleOrFile, "".__class__):
        fileHandle.close()


def getFastaDictionary(fastaFile):
    """Returns a dictionary of the first words of fasta headers to their corresponding
    fasta sequence
    """
    namesAndSequences = [(x[0].split()[0], x[1]) for x in fastaRead(open(fastaFile, 'r'))]
    names = [x[0] for x in namesAndSequences]
    assert len(names) == len(set(names))  # Check all the names are unique
    return dict(namesAndSequences)  # Hash of names to sequences

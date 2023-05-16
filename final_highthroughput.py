# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 12:14:44 2023

@author: sanga
"""

# Import Libraries
import os 
from Bio import SeqIO
from collections import defaultdict


# Set working directory
os.chdir(r'C:\Users\sanga\OneDrive\Desktop\Course Work\Spring 2023\Computational Highthroughput\Assignments')



# Load the file with RNA sequences
RNA_sequece_1 = "RNA-sequence1.fna"
RNA_sequece_2 = "RNA-sequence2.fna"

# Caluculate freq and percent of nucleotide after rtDNA
def nucleotide_freq_percent(filename, nucleotide_freq):
    """

    Parameters
    ----------
    filename : fasta file
        Takes the RNA sequences and reverse transcribe to DNA.
    nucleotide_freq : Numeric value
        Given the nucleotide frequency, the DNA sequences will be parse
        and compute frequency and percentage of the nucleotides

    Returns
    -------
    Prints sequence, freq abs and percent 
    Returns total percentages

    """
    
    sequences = []
    with open(filename, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            sequences.append(record.seq)
    print('\n' + filename , 'with', str(nucleotide_freq), "nucleotide frequency", '\n')
    # Convert RNA to DNA
    sequences_dna = [seq.back_transcribe() for seq in sequences]
    print("DNA 5'->3' sequence:")
    print(sequences_dna)
    # Calculate given nucleotide freq
    nucleotide_counts = defaultdict(int)
    
    for seq in sequences_dna:
        # Calculate given nucleotide freq
        for i in range(len(seq) - (nucleotide_freq - 1)):
            nucleotide = seq[i:i+nucleotide_freq]
            nucleotide_counts[nucleotide] += 1
    
    # Calculate total given nucleotide freq counts
    total_nucleotide_count = sum(nucleotide_counts.values())
    
    # Calculate given nucleotide freq as percent
    nucleotide_freq_percentages = {k: v/total_nucleotide_count*100 for k, v in nucleotide_counts.items()}
    
    # Print the results
    print('\n'+ str(nucleotide_freq),"Nucleotide Frequencies:")
    for nucleotide, count in nucleotide_counts.items():
        print(f"{nucleotide}: {count} (Absolute), {nucleotide_freq_percentages[nucleotide]:.2f}% (Percentage)")
    return nucleotide_freq_percentages

# The question says dinucleotide or trinucleotide. 
#I will do dinucleotide for RNA_sequece_1 and trinucleotide for RNA_sequnce_2
nucleotide_freq_percent(RNA_sequece_1, 2)
nucleotide_freq_percent(RNA_sequece_2, 2)

def diff_three_fold (seq1_percentages, seq2_percentages, nucleotide_freq, fold_diff):
    """
    Parameters
    ----------
    seq1_percentages : dict with nucleotide and percent freq
    seq2_percentages : dict with nucleotide and percent freq
    nucleotide_freq : Numeric value. Ex: 2 for dinucleotide
    fold_diff : Diff between seq1 and seq2

    Returns
    -------
    Nucleotide with difference more than 3 between seq1 and seq2

    """
    # Identify dinucleotides with a three-fold difference in percentage
    three_fold_difference = []
    for nucleotide in seq1_percentages.keys():
        if nucleotide in seq2_percentages.keys():
            percentage1 = seq1_percentages[nucleotide]
            percentage2 = seq2_percentages[nucleotide]
            if percentage1 > 0 and percentage2 > 0 and abs(percentage1 - percentage2) >= fold_diff:
                three_fold_difference.append(nucleotide)
    print('\n' + str(nucleotide_freq) , " freq nucleotides with a three-fold difference in percentage between seq1 and seq2:")
    print(three_fold_difference)

diff_three_fold(nucleotide_freq_percent(RNA_sequece_1, nucleotide_freq=2), nucleotide_freq_percent(RNA_sequece_2, nucleotide_freq=2), nucleotide_freq=2, fold_diff=3)
diff_three_fold(nucleotide_freq_percent(RNA_sequece_1, nucleotide_freq=3), nucleotide_freq_percent(RNA_sequece_2, nucleotide_freq=3), nucleotide_freq=3, fold_diff=3)

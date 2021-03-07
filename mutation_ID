# import required packages
import re
import numpy as np
import pandas as pd
from pandas import DataFrame

##### import sequences for analysis : reference and mutated

# import txt file for reference sequence
# print optional but recommended
ref_sequence_raw = open(r"C:\Users\armerv\Documents\SWBio_python_project\eGFP.txt")
ref_sequence_raw = ref_sequence_raw.read()
# print optional
#print("Reference: "+ref_sequence_raw)

# import txt file for mutated sequence
# print optional but recommended
mut_sequence_raw = open(r"C:\Users\armerv\Documents\SWBio_python_project\ZtGFP.txt")
mut_sequence_raw = mut_sequence_raw.read()
# print optional
#print("Mutated: "+mut_sequence_raw)


# find first ATG (start codon) in sequence and remove sequence prior to it
# ensures amino acid sequence in subsequent analysis and that data is suitable

# Debugging: this step was added to enable imperfect data inputs

def start_codon(seq, start):
    index = seq.find(start)
    if index !=-1 :
        return seq[index:]
    else :
        raise Exception ("Start codon not found, re-evaluate data input")


# print optional, but recommended
ref_start = start_codon(ref_sequence_raw, "ATG")
print("Reference: "+ref_start)

mut_start = start_codon(mut_sequence_raw, "ATG")
print("Mutated: "+mut_start)

# reassign raw sequence name to identified start of coding sequence.
ref_sequence_raw = ref_start
mut_sequence_raw = mut_start


##### alignment tool to align at greatest similarity
# Debugging: this should be used as a check only, to ensure data input is at start of open reading frame and hasn't identified the wrong start codon

def score(a, b):
    return sum([a[i] == b[i] for i in range(len(a))])

def best_match(ref_sequence_raw, mut_sequence_raw):
    max = 0
    max_score = 0
    for offset in range(len((ref_sequence_raw))):
        val = score(ref_sequence_raw[offset:], mut_sequence_raw)
        if val > max_score:
            max_score = val
            max = offset
    return max

alignment = best_match(ref_sequence_raw, mut_sequence_raw)
# print recommended
print("alignment score")
print(alignment)

# if alignment value is high, check data input e.g. 36 means alignment starts at base 37 (first base is 0) on the reference sequence and it is likely that the true start codon for the coding sequence has been missed 


#### Create position index with codon reference
# Debugging: this enables merging of tables after identification of changes

# number each amino acid in sequence, so not reliant on built-in index
def index(ref_sequence_raw):
    counter = 0
    Nucleotide_number = []
    Nucleotide = []
    for each_nucleotide in (ref_sequence_raw):
        Nucleotide_number.append(counter)
        Nucleotide.append(each_nucleotide)
        counter += 1
    d = ({"nucleotide_position": Nucleotide_number, "Nucleotide":Nucleotide})
    table = DataFrame(d)

    return table

Nucleotide_index_ref = index(ref_sequence_raw)
Nucleotide_index_mut = index(mut_sequence_raw)
# print optional 
#print(Nucleotide_index_ref)
#print(Nucleotide_index_mut)


# add codon position to index dataframe

def index_codon(Nucleotide_index_ref):
    
    df = pd.DataFrame(Nucleotide_index_ref)
    group_size = 3
    numbers = list(range(len(df.index) // group_size)) * group_size
    numbers.sort()
    numbers = pd.Series(numbers)
    df = pd.concat([df, numbers], ignore_index = False, axis = 1)
    df.columns = ['nucleotide_position', 'nucleotide', 'codon_position']

    groups = df.groupby('codon_position').filter(lambda x: len(x) == group_size)
    
    return groups

Index_ref = index_codon(Nucleotide_index_ref)
Index_mut = index_codon(Nucleotide_index_mut)
# print optional
#print(Index_ref)
#print(Index_mut)


##### input sequences are converted (in order) to: antisense, mRNA, codons, amino acid (3 letter) and amino acid (1 letter). Intermediary steps enable check prints if necessary and allows the user to specialise outputs to required needs.


##### convert sense strand to antisense (for mRNA)
# dictionary
sense_antisense = {
    "A":"T", "T":"A", "G":"C", "C":"G"
}

# for reference or mutated sequence
def encode_antisense(sense_to_antisense, sequence_raw):
    antisense = []

    for nucleotide in sequence_raw:
        antisense_nucleotide = sense_antisense[nucleotide]
        antisense.append(antisense_nucleotide)
        
    antisense = "".join(antisense)
    
    return antisense

# print optional
ref_antisense_seq = encode_antisense(sense_antisense, ref_sequence_raw)
#print(ref_antisense_seq)

mut_antisense_seq = encode_antisense(sense_antisense, mut_sequence_raw)
#print(mut_antisense_seq)



##### convert antisense DNA to mRNA
# dictionary
antisense_mRNA = {
    "T":"A", "A":"U", "C":"G", "G":"C"
}


# for reference sequence
def encode_mRNA(antisense_to_mRNA, antisense_seq):
    mRNA = []

    for antisense_nucleotide in antisense_seq:
        mRNA_nucleotide = antisense_mRNA[antisense_nucleotide]
        mRNA.append(mRNA_nucleotide)
    
    mRNA = "".join(mRNA)
    
    return mRNA


# print optional
ref_mRNA_seq = encode_mRNA(antisense_mRNA, ref_antisense_seq)
#print(ref_mRNA_seq)

mut_mRNA_seq = encode_mRNA(antisense_mRNA, mut_antisense_seq)
#print(mut_mRNA_seq)


##### split code into segments of 3

def encode_codon_segments(mRNA_seq):
        n = 3
        segments = [mRNA_seq[i:i+n] for i in range(0, len(mRNA_seq), n)]
        
        return segments


# print optional
ref_mRNA_segmented_seq = encode_codon_segments(ref_mRNA_seq)
#print(ref_mRNA_segmented_seq)

mut_mRNA_segmented_seq = encode_codon_segments(mut_mRNA_seq)
#print(mut_mRNA_segmented_seq)

##### convert mRNA code to 3-letter amino acid sequence (+ stop/start codons)

# dictionary
Universal_code_aa = {
    "UUU":"PHE", "UUC":"PHE", "UUA":"LEU", "UUG":"LEU", "UCU":"SER", "UCC":"SER", "UCA":"SER", "UCG":"SER",
    "UAU":"TYR", "UAC":"TYR", "UAA":"STP", "UAG":"STP", "UGU":"CYS", "UGC":"CYS", "UGA":"STP",
    "UGG":"TRP", "CUU":"LEU", "CUC":"LEU", "CUA":"LEU", "CUG":"LEU", "CCU":"PRO", "CCC":"PRO",
    "CCA":"PRO", "CCG":"PRO", "CAU":"HIS", "CAC":"HIS", "CAA":"GLN", "CAG":"GLN", "CGU":"ARG",
    "CGC":"ARG", "CGA":"ARG", "CGG":"ARG", "AUU":"ILE", "AUC":"ILE", "AUA":"ILE", "AUG":"MET",
    "ACU":"THR", "ACC":"THR", "ACA":"THR", "ACG":"THR", "AAU":"ASN", "AAC":"ASN", "AAA":"LYS",
    "AAG":"LYS", "AGU":"SER", "AGC":"SER", "AGA":"ARG", "AGG":"ARG", "GUU":"VAL", "GUC":"VAL",
    "GUA":"VAL", "GUG":"VAL", "GCU":"ALA", "GCC":"ALA", "GCA":"ALA", "GCG":"ALA", "GAU":"ASP",
    "GAC":"ASP", "GAA":"GLU", "GAG":"GLU", "GGU":"GLY", "GGC":"GLY", "GGA":"GLY", "GGG":"GLY",
}

# It is important that data input starts at a start codon for correct seperation of codons to be read as correct amino acid 
#STP is stop codon


# function that converts codons 
def encode_aa_3(Universal_code_aa, mRNA_seq):
    aa_3 = []

    for mRNA in mRNA_seq:
        amino_acid_3 = Universal_code_aa[mRNA]
        aa_3.append(amino_acid_3)
    
    aa_3 = "".join(aa_3)
    
    return aa_3

# print optional 
ref_aa_3_seq = encode_aa_3(Universal_code_aa, ref_mRNA_segmented_seq)
#print(ref_aa_3_seq)

mut_aa_3_seq = encode_aa_3(Universal_code_aa, mut_mRNA_segmented_seq)
#print(mut_aa_3_seq)

# create index dataframe for amino acid position


def index(ref_mRNA_segmented_seq):
    counter = 0
    Codon_number = []
    Codon = []
    for each_codon in ref_mRNA_segmented_seq:
        Codon_number.append(counter)
        Codon.append(each_codon)
        counter += 1
        
    d = ({"Codon":Codon, "Codon_pos_ref": Codon_number})
    table= DataFrame(d)
    
    return table   

Codon_index_ref= index(ref_mRNA_segmented_seq)
Codon_index_mut = index(mut_mRNA_segmented_seq)
#print(Codon_index_mut)
#print(Codon_index_ref)


##### segment sequences into groups of 3, reusing function from codon seperation 'encode_codon_segments'
# print optional
ref_aa_3_segmented_seq = encode_codon_segments(ref_aa_3_seq)
#print(ref_aa_3_segmented_seq)

mut_aa_3_segmented_seq = encode_codon_segments(mut_aa_3_seq)
#print(mut_aa_3_segmented_seq)


##### convert mRNA-based 3-letter amino acid sequence to a shorter, single letter amino acid sequence
# dictionary
aa_single_letter = {
    "ALA":"A", "CYS":"C", "ASP":"D", "GLU":"E", "PHE":"F", "GLY":"G", "HIS":"H", "ILE":"I", 
    "LYS":"K", "LEU":"L", "MET":"M", "ASN":"N", "PRO":"P", "GLN":"Q", "ARG":"R", "SER":"S",
    "THR":"T", "VAL":"V", "TRP":"W", "TYR":"Y", "STP":"/"
}

# function that converts 3-letter amino acid sequence to single letter amino acid sequence
def encode_aa_1(aa_single_letter, aa_3_seq):
    aa_1 = []

    for amino_acid_3 in aa_3_seq:
        amino_acid_1 = aa_single_letter[amino_acid_3]
        aa_1.append(amino_acid_1)
    
    aa_1 = "".join(aa_1)
    
    return aa_1


#print optional
ref_aa_1_seq = encode_aa_1(aa_single_letter, ref_aa_3_segmented_seq)
#print(ref_aa_1_seq)

mut_aa_1_seq = encode_aa_1(aa_single_letter, mut_aa_3_segmented_seq)
#print(mut_aa_1_seq)


##### comparison of nucleotide (DNA) sequence
# this function identifies differences between the reference and mutated sequence and outputs into a dataframe


# define variables
dna_seq_a = ref_sequence_raw
dna_seq_b = mut_sequence_raw
dna_len_a = len(dna_seq_a)    
dna_len_b = len(dna_seq_b)    


def encode_dna_sequence_compare(dna_seq_a, dna_seq_b, dna_len_a, dna_len_b):
    nucleotide_differences_a = []
    nucleotide_differences_b = []
    positions = []
    
    for nucleotide in range(0, min(dna_len_a, dna_len_b)):
        if dna_seq_a[nucleotide] != dna_seq_b[nucleotide]:
                positions.append(nucleotide)
                nucleotide_differences_a.append(dna_seq_a[nucleotide])
                nucleotide_differences_b.append(dna_seq_b[nucleotide])
    
    d = ({"nucleotide_position":positions, "ref_nucleotide":nucleotide_differences_a,
          "mut_nucleotide":nucleotide_differences_b})
    table= DataFrame(d)
    
    return table


nucleotide_seq_compare = encode_dna_sequence_compare(dna_seq_a, dna_seq_b, dna_len_a, dna_len_b)

# print optional 
#print(nucleotide_seq_compare)


#### Finding differences in amino acid sequence
# function to find differences in 3-letter amino acid 

# define variables
aa_3_seq_a = ref_aa_3_segmented_seq
aa_3_seq_b = mut_aa_3_segmented_seq
aa_3_len_a = len(aa_3_seq_a)    
aa_3_len_b = len(aa_3_seq_b) 


def encode_aa_3_sequence_compare(aa_3_seq_a, aa_3_seq_b, aa_3_len_a, aa_3_len_b):
    aa_3_differences_a = []
    aa_3_differences_b = []
    positions = []
    
    for aa_3 in range(0, min(aa_3_len_a, aa_3_len_b)):
        if aa_3_seq_a[aa_3] != aa_3_seq_b[aa_3]:
                positions.append(aa_3)
                aa_3_differences_a.append(aa_3_seq_a[aa_3])
                aa_3_differences_b.append(aa_3_seq_b[aa_3])
    
    d = ({"codon_position":positions, "ref_aa":aa_3_differences_a, "mut_aa":aa_3_differences_b})
    table= DataFrame(d)
    
    return table


aa_3_seq_compare = encode_aa_3_sequence_compare(ref_aa_3_segmented_seq, mut_aa_3_segmented_seq, aa_3_len_a, aa_3_len_b)
# print optional
#print(aa_3_seq_compare)


#### codon usage 

codon_seq_a = ref_mRNA_segmented_seq
codon_seq_b = mut_mRNA_segmented_seq
codon_len_a = len(codon_seq_a)    
codon_len_b = len(codon_seq_b) 


def encode_codon_sequence_compare(codon_seq_a, codon_seq_b, codon_len_a, codon_len_b):
    codon_differences_a = []
    codon_differences_b = []
    positions = []
    
    for codon in range(0, min(codon_len_a, codon_len_b)):
        if codon_seq_a[codon] != codon_seq_b[codon]:
                positions.append(codon)
                codon_differences_a.append(codon_seq_a[codon])
                codon_differences_b.append(codon_seq_b[codon])
    
    d = ({"codon_position":positions, "ref_codon":codon_differences_a, "mut_codon":codon_differences_b})
    table= DataFrame(d)
    # codon position 0 is first i.e. ATG / start codon
    return table


codon_seq_compare = encode_codon_sequence_compare(ref_mRNA_segmented_seq, mut_mRNA_segmented_seq, codon_len_a, codon_len_b)
# print optional
#print(codon_seq_compare)


##### integrate all results into single table 
# merge data tables
# delete unnecessary column
del Index_ref['nucleotide']

# merge of tables done step-wise to ensure correct alignment
# merge index and nucleotide differences
merged_results = pd.merge(left=Index_ref, right= nucleotide_seq_compare, left_on = "nucleotide_position", right_on = "nucleotide_position")
#print(merged_results)

# add codon differences
Results_1 = pd.merge(left = merged_results, right = codon_seq_compare, left_on = "codon_position", right_on = "codon_position")

# add amino acid differences
Results_2 = pd.merge(left = Results_1, right = aa_3_seq_compare, how="left", left_on = "codon_position", right_on = "codon_position")
# reorder columns
Results = Results_2[["nucleotide_position", "ref_nucleotide", "mut_nucleotide", "codon_position", "ref_codon", "mut_codon", "ref_aa", "mut_aa"]]
# replace NAN with blank space in dataframe
Results_final = Results.replace(np.nan, "", regex= True)
# rename columns (to fit onto single line) -optional
Results_final = Results_final.rename(columns={"nucleotide_position": "n_pos", "ref_nucleotide": "ref_n", "mut_nucleotide": "mut_n", "codon_position":"c_pos"})

#### final outputs

# final results table 
print("Identified mutations")
print(Results_final)


##### summary table 

# Number of mutations count
num_mutations = nucleotide_seq_compare.count()
num_mutations = DataFrame(num_mutations)
num_mutations = num_mutations.T


# number of amino acid changes count
num_aa_changes = aa_3_seq_compare.count()
num_aa_changes = DataFrame(num_aa_changes)
num_aa_changes = num_aa_changes.T


# delete unnecessary columns
del num_mutations["ref_nucleotide"]
del num_mutations["mut_nucleotide"]
del num_aa_changes["ref_aa"]
del num_aa_changes["mut_aa"]

# merge tables
Summary_table = pd.concat([num_mutations, num_aa_changes], axis =1).reindex(num_mutations.index)
# rename colums
Summary_table = Summary_table.rename(columns={"nucleotide_position": "nucleotide_changes", "codon_position":"amino_acid_changes"})

#### summary of mutations and amino acid changes

# print summary table 
print("Summary Table")
print(Summary_table)


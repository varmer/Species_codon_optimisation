# Victoria Armer
# for nucleotide optimisation


# This code calculates codon frequency usage from predicted coding sequences of the Fusarium graminearum genome and converts protein sequences into codon-optimzed amino acid sequences for increased transcription and translation by the host organism. This is useful for transgenic studies. Examples here are GFP and mCherry for use in microscopy.



# Set up: download cds of host organism from ensembl. Fusarium graminearum genome can be downloaded here: ftp://ftp.ensemblgenomes.org/pub/fungi/release-49/fasta/fusarium_graminearum/cds/ (Feb, 2021)

# import required packages

import pandas as pd
from pandas import DataFrame
import re
import numpy as np


#### PART 1 
# Determining species codon usage 



# import coding sequences of fusarium 
### change so that this is a file path selector before uploading to github
Species_cds = open(r"C:\Users\armerv\Desktop\BLAST\Genome_cds_fastas\Fusarium_graminearum_cds.fasta")

#Species_cds_test = ['ATGTATTCGTGCCTAGCATGACTACGATGACATGACAATGCGCA']

#### ignore cds descriptions starting with '>'

# this function removes the coding sequence description line starting '>' to return a list of just the coding sequences
# species is the file path, which we have already defined as 'Species_cds' for 'Species coding sequences'

def get_coding_seqs(cds):
    seqs = []
    for line in cds:
        if line.startswith('>') is False:
            seqs.append(line.rstrip())
    return seqs


# convert DNA to mRNA in truncated coding sequences

def DNA_to_mRNA(seqs):
    sequences = []
    for cds in seqs:
        sequences.append(cds.upper().replace('T', 'U'))    
    return sequences


# split coding seqs into codons (all start at start codon already, but remove any that aren't factor of 3)

def codon_segmentation(sequences):
        n = 3
        segments = []
        for seq in sequences:
            for i in range(0, len(seq), n):
                segments.append(seq[i:i+n])
        return segments

# create a dictionary of amino acid codons
##break up 6-codon family into 2 and 4 fold (Ser (S), L(Leu), R (Arg))
c_to_aa = {
    "UUU":"Phe", "UUC":"Phe",         
    "UCU":"Ser4", "UCC":"Ser4", "UCA":"Ser4", "UCG":"Ser4",
    "AGU":"Ser2", "AGC":"Ser2",
    "CUU":"Leu4", "CUC":"Leu4", "CUA":"Leu4", "CUG":"Leu4",
    "UUA":"Leu2", "UUG":"Leu2",
    
    "UAU":"Tyr", "UAC":"Tyr", "UAA":"STOP", "UAG":"STOP",
    "UGU":"Cys", "UGC":"Cys", "UGA":"STOP", "UGG":"Trp",
    "CGU":"Arg4", "CGC":"Arg4", "CGA":"Arg4", "CGG":"Arg4",
    "AGA":"Arg2", "AGG":"Arg2",
    "CCU":"Pro", "CCC":"Pro", "CCA":"Pro", "CCG":"Pro",
    "CAU":"His", "CAC":"His", "CAA":"Gln", "CAG":"Gln",
    
    "AUU":"Ile", "AUC":"Ile", "AUA":"Ile", "AUG":"Met",
    "ACU":"Thr", "ACC":"Thr", "ACA":"Thr", "ACG":"Thr",
    "AAU":"Asn", "AAC":"Asn", "AAA":"Lys", "AAG":"Lys",
   
    "GUU":"Val", "GUC":"Val", "GUA":"Val", "GUG":"Val",
    "GCU":"Ala", "GCC":"Ala", "GCA":"Ala", "GCG":"Ala",
    "GAU":"Asp", "GAC":"Asp", "GAA":"Glu", "GAG":"Glu",
    "GGU":"Gly", "GGC":"Gly", "GGA":"Gly", "GGG":"Gly"}


#### Retrieving codon frequencies
# seqs = list of truncadted mRNA coding sequences
# returns dataframe of absolute codon frequencies


def counter(segments):
    codon_count = dict()
    
    for codon in list(c_to_aa.keys()): # use dictionary for 64 codons
        codon_count[codon]=0  # empty vector to accumulate codon count

    for c in list(codon_count.keys()):
        codon_count[c]+= segments.count(c)

    df_codon_count=pd.DataFrame(list(codon_count.items()) ) #output as dataframe
    df_codon_count.columns=['Codon', 'Obs_Freq'] # column names
    df_codon_count['Amino_Acid'] = [c_to_aa[codon] for codon in df_codon_count['Codon'].values] ## add amino acid column

    return df_codon_count

#### caluculate codon usage from codon counts

def codon_usage(df_codon_count): #pass in codon count from previous
    aa_groups = df_codon_count.groupby('Amino_Acid') # by amino acid column in dataframe
    aa =  df_codon_count['Amino_Acid'].unique()  #make a list of all amino acids to iterate over
    df_list = []
    for a in aa:
        d=aa_groups.get_group(a)
        d['Relative_codon_usage'] = (d['Obs_Freq'].values)/(d['Obs_Freq'].mean()) #obs/expected freq 
        d['Relative_adaptive_weights']= (d['Relative_codon_usage'].values)/(d['Relative_codon_usage'].max())
        d['optimal'] = [True if relative_usage ==d['Relative_codon_usage'].max() else False for relative_usage in d['Relative_codon_usage'].values] 
        #marks optimal codon for each aa
        df_list.append(d) # creates lists

    return pd.concat(df_list) # joins lists



#### Run programmes sequentially, putting output of previous as input
truncation = get_coding_seqs(Species_cds)
mRNA_conversion = DNA_to_mRNA(truncation)
segmentation = codon_segmentation(mRNA_conversion)
codon_counter = counter(segmentation)
codon_bias = codon_usage(codon_counter)
pd.set_option('display.max_columns', None) # allows all columns to be displayed
print(codon_bias)


##saves final table as csv (if required)
#codon_bias.to_csv('Species_codon_usage.csv')



#### PART 2
# optimizing protein sequence according to species codon usage bias


# import protein sequence for desired optimization. Here eGFP and mCherry are imported.
protein_cds = open('file path here')

# run previous function to parse sequence, convert to mRNA and segment
protein_truncation = get_coding_seqs(protein_cds)
protein_mRNA_conversion = DNA_to_mRNA(protein_truncation)
protein_segmentation = codon_segmentation(protein_mRNA_conversion)


# convert sequence to amino acid sequence using dictionary

# function that converts codons - NOT YET TESTED
def amino_acid_conversion(c_to_aa, seg_protein_seq):
    aa = []

    for protein_codons in seg_protein_seq:
        amino_acid = c_to_aa[protein_codons]
        aa.append(amino_acid)
    
    aa = "".join(aa)
    
    return aa

protein_aa = amino_acid_conversion(protein_segmentation)
# check
print(protein_aa)


#### need to pull out list of top codon for each amino acid - already been marked in csv 'Species_codon_usage'
# delete all calculation columns are rows where optimal = FALSE 

df = codon_bias # create new table
del df['Obs_Freq']
del df['Relative_codon_usage']
del df['Relative_adaptive_weights']
   
#check     
print(df)

# filter rows for optimal codon
df = df[df.optimal != False]

#remove optimal column as no longer required
del df['optimal']

# create dictionary converting amino acids to codons 
Top_codons = dict(df.values)
print(Top_codons)


# function to convert aa back to optimal codons

def codon_conversion(Top_codons, protein_aa):
    new_codons = []

    for amino_acid in protein_aa:
        new_codon = Top_codons[protein_aa]
        new_codons.append(new_codon)
    
    new_codons= "".join(new_codons)
    
    return new_codons


# run function
protein_new_cds = codon_conversion(protein_aa)
#check
print(protein_new_cds)



### output sequence, saved as fasta file


    

### sequence comparison i.e. print differences
    
  
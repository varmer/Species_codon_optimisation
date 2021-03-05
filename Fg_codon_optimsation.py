# Victoria Armer
# for nucleotide optimisation


# This code calculates codon frequency usage from predicted coding sequences of the Fusarium graminearum genome and converts protein sequences into codon-optimzed amino acid sequences for increased transcription and translation by the host organism. This is useful for transgenic studies. Examples here are GFP and mCherry for use in microscopy.



# Set up: download cds of host organism from ensembl. Fusarium graminearum genome can be downloaded here: ftp://ftp.ensemblgenomes.org/pub/fungi/release-49/fasta/fusarium_graminearum/cds/ (Feb, 2021)

# import required packages

import pandas as pd
from pandas import DataFrame
import re
import numpy as np


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


# truncation check
trunc_check = get_coding_seqs(Species_cds)
#print(trunc_check)


# convert DNA to mRNA in truncated coding sequences

def DNA_to_mRNA(seqs):
    sequences = []
    for cds in seqs:
        sequences.append(cds.upper().replace('T', 'U'))    
    return sequences

# mRNA check

mRNA_check = DNA_to_mRNA(trunc_check)
#print(mRNA_check)



# split coding seqs into codons (all start at start codon already, but remove any that aren't factor of 3)

def encode_codon_segments(mRNA_seq):
        n = 3
        segments = []
        for seq in mRNA_seq:
            for i in range(0, len(seq), n):
                segments.append(seq[i:i+n])
        return segments

# segmentation check
seg_check = encode_codon_segments(mRNA_check)
#print(seg_check)



# create a dictionary of amino acid codons
##break up 6-codon family into 2 and 4 fold (Ser (S), L(Leu), R (Arg))
codon_to_aa = {
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


def counter(codons):
    codon_count = dict()
    
    for codon in list(codon_to_aa.keys()): # use dictionary for 64 codons
        codon_count[codon]=0  # empty vector to accumulate codon count

    for c in list(codon_count.keys()):
        codon_count[c]+= codons.count(c)

    df_codcount=pd.DataFrame(list(codon_count.items()) ) #output as dataframe
    df_codcount.columns=['Codon', 'Obs_Freq'] # column names
    df_codcount['Amino_Acid'] = [codon_to_aa[codon] for codon in df_codcount['Codon'].values] ## add amino acid column

    return df_codcount

Count_check = counter(seg_check)
print(Count_check)


#### caluculate codon usage from codon counts

def calculate_codon_usage(df_codcount): #pass in codon count from previous
    aa_groups = df_codcount.groupby('Amino_Acid') # by amino acid column in dataframe
    aa =  df_codcount['Amino_Acid'].unique()  #make a list of all amino acids to iterate over
    df_list = []
    for a in aa:
        d=aa_groups.get_group(a)
        d['Relative_codon_usage'] = (d['Obs_Freq'].values)/(d['Obs_Freq'].mean()) #obs/expected freq 
        d['Relative_adaptive_weights']= (d['Relative_codon_usage'].values)/(d['Relative_codon_usage'].max())
        d['optimal'] = [True if relative_usage ==d['Relative_codon_usage'].max() else False for relative_usage in d['Relative_codon_usage'].values] 
        #marks optimal codon for each aa
        df_list.append(d) # creates lists

    return pd.concat(df_list) # joins lists

#check
Usage = calculate_codon_usage(Count_check)
pd.set_option('display.max_columns', None)
print(Usage)



#### Run programmes sequentially, putting output of previous as input
coding_seqs = get_coding_seqs(Species_cds)
codon_freq_check = acquire_codon_freq(coding_seqs)
codon_usage = calculate_codon_usage(codon_freq_check)

##saves final table as csv

codon_usage.to_csv('Species_codon_usage.csv')



#### need to pull out list of top codon for each amino acid - already been marked in csv 'Species_codon_usage'










# import protein sequence for desired optimization. Here eGFP and mCherry are imported.
#eGFP = open
#mCherry = open





### parse cds : remove '>'




### identify codons in sequence




### replace function to put in top codon




### output sequence, saved as fasta file


    

### sequence comparison i.e. print differences
    
  
# Written by Victoria Armer Feb 2021
# A script for determining species codon bias and protein sequence codon optimzation 

# Set up: download cds of host organism from ensembl. Fusarium graminearum genome can be downloaded here: ftp://ftp.ensemblgenomes.org/pub/fungi/release-49/fasta/fusarium_graminearum/cds/ (Feb, 2021)

# very quick  ~4 seconds for Fusarium graminearum genome

# import required packages
import pandas as pd

#### PART 1 
# Determining species codon bias

# import coding sequences of species (e.g. Fusarium graminearum)
### change so that this is a file path selector before uploading to github
Species_cds = open(r"C:\Users\armerv\Desktop\BLAST\Genome_cds_fastas\Fusarium_graminearum_cds.fasta")

# remove first line from fasta sequences with descriptions

def get_coding_seqs(cds):
    code = []
    for line in cds:
        if line.startswith('>') is False:
            code.append(line.rstrip())
    return code


# convert DNA to mRNA in truncated coding sequences

def DNA_to_mRNA(DNA):
    mRNA = []
    for code in DNA:
        mRNA.append(code.upper().replace('T', 'U'))    
    return mRNA


# split coding seqs into codons (all start at start codon already, but remove any that aren't factor of 3)

def codon_segmentation(mRNA):
        n = 3
        segments = []
        for seq in mRNA:
            for i in range(0, len(seq), n):
                segments.append(seq[i:i+n])
        return segments

# create a dictionary of amino acid codons
c_to_aa = {
    "UUU":"Phe", "UUC":"Phe",         
    "UCU":"Ser", "UCC":"Ser", "UCA":"Ser", "UCG":"Ser",
    "AGU":"Ser", "AGC":"Ser",
    "CUU":"Leu", "CUC":"Leu", "CUA":"Leu", "CUG":"Leu",
    "UUA":"Leu", "UUG":"Leu",
    "UAU":"Tyr", "UAC":"Tyr", "UAA":"STOP", "UAG":"STOP",
    "UGU":"Cys", "UGC":"Cys", "UGA":"STOP", "UGG":"Trp",
    "CGU":"Arg", "CGC":"Arg", "CGA":"Arg", "CGG":"Arg",
    "AGA":"Arg", "AGG":"Arg",
    "CCU":"Pro", "CCC":"Pro", "CCA":"Pro", "CCG":"Pro",
    "CAU":"His", "CAC":"His", "CAA":"Gln", "CAG":"Gln",
    "AUU":"Ile", "AUC":"Ile", "AUA":"Ile", "AUG":"Met",
    "ACU":"Thr", "ACC":"Thr", "ACA":"Thr", "ACG":"Thr",
    "AAU":"Asn", "AAC":"Asn", "AAA":"Lys", "AAG":"Lys",
    "GUU":"Val", "GUC":"Val", "GUA":"Val", "GUG":"Val",
    "GCU":"Ala", "GCC":"Ala", "GCA":"Ala", "GCG":"Ala",
    "GAU":"Asp", "GAC":"Asp", "GAA":"Glu", "GAG":"Glu",
    "GGU":"Gly", "GGC":"Gly", "GGA":"Gly", "GGG":"Gly"}


# codon incidences in genome

def counter(segments):
    codon_count = dict() # create dictionary
    
    for codon in list(c_to_aa.keys()): # use dictionary for codons
        codon_count[codon]=0  # empty vector to accumulate codon count

    for c in list(codon_count.keys()):
        codon_count[c]+= segments.count(c) #counter

    df_codon_count=pd.DataFrame(list(codon_count.items()) ) #output as dataframe
    df_codon_count.columns=['Codon', 'Count'] # column names
    df_codon_count['Amino_Acid'] = [c_to_aa[codon] for codon in df_codon_count['Codon'].values] ## add amino acid column

    return df_codon_count

# caluculate relative codon usage from codon incidence

def codon_usage(df_codon_count): #pass in codon count from previous
    aa_groups = df_codon_count.groupby('Amino_Acid') # by amino acid column in dataframe
    aa =  df_codon_count['Amino_Acid'].unique()  #make a list of all amino acids to iterate over
    df_list = []
    for a in aa:
        d=aa_groups.get_group(a)
        d['Relative_codon_usage'] = (d['Count'].values)/(d['Count'].mean()) #obs/expected freq 
        d['Relative_weighting']= (d['Relative_codon_usage'].values)/(d['Relative_codon_usage'].max())
        d['Choice_codon'] = [True if relative_usage ==d['Relative_codon_usage'].max() else False for relative_usage in d['Relative_codon_usage'].values] 
        #marks optimal codon for each aa
        df_list.append(d) # creates lists

    return pd.concat(df_list) # joins lists



#### Run programmes sequentially

truncation = get_coding_seqs(Species_cds)
mRNA_conversion = DNA_to_mRNA(truncation)
segmentation = codon_segmentation(mRNA_conversion)
codon_counter = counter(segmentation)
codon_bias = codon_usage(codon_counter)
pd.set_option('display.max_columns', None) # allows all columns to be displayed
print(codon_bias) # print results


##saves final table as csv (if required)
#codon_bias.to_csv('Species_codon_usage.csv')



#### PART 2
# optimizing protein sequence according to species codon usage bias


# import protein sequence for desired optimization. Here eGFP is imported.
protein_cds =  open(r"C:\Users\armerv\Desktop\BLAST\Protein_seqs\eGFP.fasta")
#protein_cds = protein_cds.read()

# run previous function to parse sequence
protein_truncation = get_coding_seqs(protein_cds)

# split string into single amino acid sequence
def single_aa_segmentation(protein_truncation):
        n = 1
        segments = []
        for seq in protein_truncation:
            for i in range(0, len(seq), n):
                segments.append(seq[i:i+n])
        return segments

# run
protein_segmentation = single_aa_segmentation(protein_truncation)


# convert single letter amini acid sequence to 3 letter aa sequence using dictionary
aa_single_letter = {
    "A":"Ala", "C":"Cys", "D": "Asp", "E":"Glu", "F":"Phe", "G":"Gly", "H":"His", "I":"Ile", 
    "K":"Lys", "L":"Leu", "M":"Met", "N":"Asn", "P":"Pro", "Q":"Gln", "R":"Arg", "S":"Ser",
    "T":"Thr", "V":"Val", "W":"Trp", "Y":"Tyr"
}


# Convert codons
def amino_acid_conversion(aa_single_letter, protein_segmentation):
    aa = []

    for single_aa in protein_segmentation:
        amino_acid = aa_single_letter[single_aa]
        aa.append(amino_acid)
    
    
    return aa

protein_aa = amino_acid_conversion(aa_single_letter, protein_segmentation)

# reduce table to just include top codons for each amino acid

df = codon_bias # create new table
# delete unnecessary columns
del df['Count']
del df['Relative_codon_usage']
del df['Relative_weighting']

# filter rows for optimal codon
df = df[df.Choice_codon != False]

#remove optimal column as no longer required
del df['Choice_codon']


# reorder columns so dictionary is right way
df = df[["Amino_Acid", "Codon"]]
# create dictionary converting amino acids to codons 
Top_codons = dict(df.values)

# convert aa to optimal codons for species
def codon_conversion(Top_codons, protein_aa):
    new_codons = []

    for amino_acid in protein_aa:
        new_codon = Top_codons[amino_acid]
        new_codons.append(new_codon)
    
    new_codons= "".join(new_codons)
    
    return new_codons


# run function
protein_new_cds = codon_conversion(Top_codons, protein_aa)

# segment mRNA sequence to replace U with T for DNA sequence
protein_new_cds_1 = single_aa_segmentation(protein_new_cds)

# convert back to DNA
def mRNA_to_DNA(protein_new_cds_1):
    new_DNA = []
    
    for cds in protein_new_cds_1:
        new_DNA.append(cds.upper().replace('U', 'T'))
        
    new_DNA = "".join(new_DNA)
    
    return new_DNA


optimized_protein = mRNA_to_DNA(protein_new_cds_1)

# check
print(optimized_protein)



    
  
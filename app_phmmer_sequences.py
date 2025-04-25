import os

from funct.phmmer import phmmer_filter_keywords, phmmer_filter_E
from funct.fetch_sequences import retrieve_phmmer_sequences, retrieve_uniprot_sequences
from funct.fasta_handling import add_fasta2fasta

import pandas as pd

# EDIT
file_prefix = "phmmer_rp15" # specify input CSV with phmmer hits
df_phmmer = pd.read_csv(f"data/phmmer_results/{file_prefix}.csv", delimiter=",", header=0) # (when Uniprot database, delimiter must be ';'


# EDIT
# choose categories and their keywords # EDIT
keywords_dict = {
            "key1": ["Dehydrogenase"],
            "key2": ["Reductase"],
            "key3": ["Oxidoreductase", "Oxido-reductase", "Oxido-Reductase"]
        }

key_split_dict = phmmer_filter_keywords(df_phmmer, f"data/phmmer_results/keyword_filter/{file_prefix}_", keywords_dict=keywords_dict)

# continue with subset of phmmer hits
subset1 = "key1" # EDIT
df_dld = key_split_dict[subset1]

# filter subset for e-values
e_filter_dict = phmmer_filter_E([df_dld], [1e-250, 1e-200, 1e-100, 1e-75, 1e-50]) # apply E-value filters # EDIT
df_dld_e = e_filter_dict[1e-50] # choose E-value result # EDIT
e = "E-50" # EDIT

# BE AWARE WHEN SETTING TO "False" !!!!
sequences_retrieved = False # RETRIEVE SEQUENCES ONLY ONCE!
# All sequences have already been retrieved, IT TAKES A LOT OF TIME!

save_fasta = f"data/phmmer_results/filtered/{file_prefix}_{subset1}_{e}.fasta"

if not sequences_retrieved:
    retrieve_phmmer_sequences(df_dld_e, save_fasta)
    #retrieve_uniprot_sequences(df_dld_e, save_fasta)

# Add query sequence to phmmer results fasta
complete_fasta = f"{save_fasta.split('.')[0]}_query.fasta"
file = complete_fasta.split("/")[-1]
aln_fasta = f"data/clustalo/{file.split('.')[0]}_aln.fasta"
add_fasta2fasta("data/query.fasta", save_fasta, complete_fasta)

# MSA:
cmd = f'C:/Users/c-fih/anaconda3/envs/clustalo/clustalo-122-win64/clustal-omega-1.2.2-win64/clustalo -i {complete_fasta} -o {aln_fasta}'
os.system(cmd)

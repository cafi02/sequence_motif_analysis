from funct.phmmer import phmmer_random, phmmer_filter_keywords
from funct.fetch_sequences import retrieve_phmmer_sequences
from funct.fasta_handling import add_fasta2fasta, fasta_to_biotite_seq
from funct.dendrogram import get_hit_description, dendrogram_from_aln

import biotite.application.clustalo as clustalo
import random
import pandas as pd


df_phmmer = pd.read_csv("data/phmmer_results/phmmer_rp15.csv", header=0) # EDIT
output_prefix = "data/phmmer_results/random/phmmer_rp15" # EDIT

new_random_file = True # EDIT
new_rand5_file = True # EDIT
new_randDLD_file = True # EDIT



db = output_prefix.split("_")[-1]
random_state = random.randint(0, 10000)

if new_random_file:
    # Select 20 random hits
    df_random = phmmer_random([df_phmmer], output_file=f"{output_prefix}_random.csv",
                              num_rows=20, seed=random_state)
    print(df_random.shape)
    retrieve_phmmer_sequences(df_random, f"{output_prefix}_random.fasta")
    add_fasta2fasta("data/RdE3.fasta", f"{output_prefix}_random.fasta", f"{output_prefix}_random_RdE3.fasta")

if new_rand5_file:
    # Select 5 random hits of each keyword-subset and concatenate
    key_split_dict = phmmer_filter_keywords(df_phmmer, f"data/phmmer_results/keyword_filter/phmmer_{db}_", keywords_dict=None) # default dict is defined in function for DLD
    df_rand5 = phmmer_random(list(key_split_dict.values()), output_file=f"{output_prefix}_rand5.csv",
                              num_rows=5, seed=random_state)
    print(df_rand5.shape)
    retrieve_phmmer_sequences(df_rand5, f"{output_prefix}_rand5.fasta")
    add_fasta2fasta("data/RdE3.fasta", f"{output_prefix}_rand5.fasta", f"{output_prefix}_rand5_RdE3.fasta")

if new_randDLD_file:
    # Select 20 random hits of the DLD-subset
    key_split_dict = phmmer_filter_keywords(df_phmmer, "data/phmmer_results/keyword_filter/phmmer_rp15_", keywords_dict=None) # default dict is defined in function for DLD
    df_randDLD = phmmer_random([key_split_dict["DLD"]], output_file=f"{output_prefix}_randDLD.csv",
                              num_rows=20, seed=random_state)
    print(df_randDLD.shape)
    retrieve_phmmer_sequences(df_randDLD, f"{output_prefix}_randDLD.fasta")
    add_fasta2fasta("data/RdE3.fasta", f"{output_prefix}_randDLD.fasta", f"{output_prefix}_randDLD_RdE3.fasta")


# Create Dendrogram from random subset of pHMMER hits
random_dict, random_seqs, random_ids, random_starts = fasta_to_biotite_seq(f"{output_prefix}_random_RdE3.fasta",
                                                                           full_headers=True)
random_ids = get_hit_description(random_ids)
random_aln = clustalo.ClustalOmegaApp.align(random_seqs,
                                            bin_path="C:/Users/c-fih/anaconda3/envs/clustalo/clustalo-122-win64/clustal-omega-1.2.2-win64/clustalo.exe")
dendrogram_from_aln(random_aln, random_ids, title="Sequence deviation of a random subset from pHMMER hits",
                    png_out=f"data/dendrograms/{db}_dendrogram_random.png")

# Create Dendrogram from 5 random hits of each keyword-subset of pHMMER hits
rand5_dict, rand5_seqs, rand5_ids, rand5_starts = fasta_to_biotite_seq(f"{output_prefix}_rand5_RdE3.fasta",
                                                                       full_headers=True)
rand5_ids = get_hit_description(rand5_ids)
rand5_aln = clustalo.ClustalOmegaApp.align(rand5_seqs,
                                            bin_path="C:/Users/c-fih/anaconda3/envs/clustalo/clustalo-122-win64/clustal-omega-1.2.2-win64/clustalo.exe")
dendrogram_from_aln(rand5_aln, rand5_ids, title="Sequence deviation of 5 random hits of each keyword-subset from pHMMER hits",
                    png_out=f"data/dendrograms/{db}_dendrogram_rand5.png")

# Create Dendrogram from 20 random hits of the DLD-subset
randDLD_dict, randDLD_seqs, randDLD_ids, randDLD_starts = fasta_to_biotite_seq(f"{output_prefix}_randDLD_RdE3.fasta",
                                                                               full_headers=True)
randDLD_ids = get_hit_description(randDLD_ids)
randDLD_aln = clustalo.ClustalOmegaApp.align(randDLD_seqs,
                                             bin_path="C:/Users/c-fih/anaconda3/envs/clustalo/clustalo-122-win64/clustal-omega-1.2.2-win64/clustalo.exe")
dendrogram_from_aln(randDLD_aln, randDLD_ids, title="Sequence deviation of 20 random hits of the DLD-subset from pHMMER hits",
                    png_out=f"data/dendrograms/{db}_dendrogram_randDLD.png")
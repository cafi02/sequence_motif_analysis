from funct.taxonomy import extract_organisms_from_fasta, get_taxonomic_distribution

# EDIT
file = "phmmer_rp15_key1_E-50_query"  # specify FASTA file, that you want to get the taxonomic distribution from

organisms = extract_organisms_from_fasta(f"data/phmmer_results/filtered/{file}.fasta")
phyla = get_taxonomic_distribution(organisms, f"data/tax_dist/{file}_dist", "phylum")
kingdom = get_taxonomic_distribution(organisms, f"data/tax_dist/{file}_dist", "kingdom")
domain = get_taxonomic_distribution(organisms, f"data/tax_dist/{file}_dist", "domain")


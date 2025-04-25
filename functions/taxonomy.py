from ete3 import NCBITaxa
import matplotlib.pyplot as plt
from collections import defaultdict

from funct.dendrogram import abbreviate


def extract_organisms_from_fasta(file_path, abb=False):

    """
    Extracts organism names from a FASTA file.

    Parameters:
        file_path (str): Path to the FASTA file.
        abb (bool): If True, returns abbreviated organism names.

    Returns:
        list: A list of extracted organism names, abbreviated if specified.
    """
    organisms = []
    with open(file_path, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):  # Header line in FASTA format
                # Split the line at the pipe "|" and extract the organism info
                header = line
                organism = header.split("OS=")[1]
                organism = organism.split("OX=")[0]
                organism = organism.strip()
                #organism = get_organism_from_fasta_header(header)
                #print(organism)
                organisms.append(organism)

    if abb:
        organisms = [abbreviate(species) for species in organisms]

    return organisms






def get_taxonomic_distribution(organisms, output_file, rank="phylum"):
    """
    Retrieves the taxonomic rank (e.g., phylum, class, order) of given organisms using the NCBI taxonomy database,
    plots the distribution as a bar chart, saves it as a PNG file, and returns a dictionary mapping taxonomic rank to organisms.

    Args:
    organisms (list): List of organism names (strings).
    output_file (str): Base name for output files (without extension).
    rank (str): Taxonomic rank to retrieve (default: "phylum").

    Returns:
    dict: Dictionary with the specified rank as keys and lists of organisms as values.
    """

    # Initialize NCBI Taxonomy database
    ncbi = NCBITaxa()

    taxonomic_dict = defaultdict(list)

    for organism in organisms:
        try:
            taxid = ncbi.get_name_translator([organism])
            if not taxid:
                print(f"Warning: Could not find taxonomy ID for {organism}")
                continue

            taxid = list(taxid.values())[0][0]  # Extract first taxonomic ID
            lineage = ncbi.get_lineage(taxid)  # Get lineage (list of taxonomic IDs)
            ranks = ncbi.get_rank(lineage)  # Get rank mapping {taxid: rank}

            # Find the taxid corresponding to the specified rank
            rank_taxid = next((tid for tid, r in ranks.items() if r == rank), None)

            if rank_taxid:
                rank_name = ncbi.get_taxid_translator([rank_taxid])[rank_taxid]
                taxonomic_dict[rank_name].append(organism)
            else:
                print(f"Warning: Could not determine {rank} for {organism}")

        except Exception as e:
            print(f"Error processing {organism}: {e}")

    # Plot Taxonomic Rank Distribution
    taxonomic_groups = list(taxonomic_dict.keys())
    counts = [len(taxonomic_dict[group]) for group in taxonomic_groups]

    size = len(taxonomic_groups) * 0.3 + 1



    plt.figure(figsize=(size, 4))
    plt.bar(taxonomic_groups, counts, edgecolor='black', color='skyblue')
    plt.xlabel(rank.capitalize(), fontweight="bold", fontsize=11)  # Capitalize for display
    plt.ylabel("Number of Organisms", fontweight="bold", fontsize=11)
    plt.title(f"{rank.capitalize()} Distribution")

    plt.xticks(rotation=45, ha="right", fontsize=9)

    plt.tight_layout()

    # Save plot
    plot_file = f"{output_file}_{rank}.png"
    plt.savefig(plot_file, dpi=300)
    print(f"Plot saved as {plot_file}")

    # Save dictionary as text file
    dict_file = f"{output_file}_{rank}.txt"
    for key, value in taxonomic_dict.items():
        print(f"{key}: {len(value)}")
    with open(dict_file, "w") as f:
        for taxon, orgs in taxonomic_dict.items():
            f.write(f"{taxon}: {', '.join(orgs)}\n")

    print(f"{rank.capitalize()} dictionary saved as {dict_file}")

    return dict(taxonomic_dict)

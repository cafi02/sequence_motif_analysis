import requests
import time
from tqdm import tqdm


def retrieve_phmmer_sequences(df, output_file):
    """
    Retrieves sequences from UniProt's online database based on IDs in a DataFrame,
    and displays a progress bar.

    Parameters:
        df (pandas.DataFrame): DataFrame containing a "Target Name" column with UniProt IDs.
        output_file (str): Path to save the retrieved sequences in FASTA format.

    Returns:
        None
    """

    base_url = "https://rest.uniprot.org/uniprotkb/"
    target_ids = set(df["Target Name"].astype(str))  # Ensure IDs are strings

    with open(output_file, "w") as fasta_out, tqdm(total=len(target_ids), desc="Fetching Sequences") as pbar:
        for uniprot_id in target_ids:
            url = f"{base_url}{uniprot_id}.fasta"

            try:
                response = requests.get(url, timeout=10)

                if response.status_code == 200:
                    fasta_out.write(response.text)
                else:
                    print(f"\n❌ Not found: {uniprot_id} (HTTP {response.status_code})")

            except requests.RequestException as e:
                print(f"\n⚠️ Error retrieving {uniprot_id}: {e}")

            pbar.update(1)  # Update the progress bar
            time.sleep(0.5)  # Prevent overloading UniProt's servers

    print(f"✅ Finished retrieving sequences. Saved to {output_file}")


import requests
import pandas as pd
import time
from tqdm import tqdm
import re  # Import regex module


def retrieve_uniprot_sequences(df, output_file):
    """
    Retrieves sequences from UniProt's online database based on IDs in a DataFrame,
    and displays a progress bar.

    Parameters:
        df (pandas.DataFrame): DataFrame containing a "Target Name" column with UniProt IDs.
        output_file (str): Path to save the retrieved sequences in FASTA format.

    Returns:
        None
    """

    base_url = "https://rest.uniprot.org/uniprotkb/"

    def extract_uniprot_id(identifier):
        """Extracts the UniProt accession number from an ID like 'sp|O08749|DLDH_MOUSE'."""
        match = re.search(r"\|([A-Z0-9]+)\|", identifier)  # Find text between '|' characters
        return match.group(1) if match else identifier  # Return found accession or original text

    target_ids = set(df["Target Name"].astype(str).apply(extract_uniprot_id))  # Extract valid IDs

    with open(output_file, "w") as fasta_out, tqdm(total=len(target_ids), desc="Fetching Sequences") as pbar:
        for uniprot_id in target_ids:
            url = f"{base_url}{uniprot_id}.fasta"

            try:
                response = requests.get(url, timeout=10)

                if response.status_code == 200:
                    fasta_out.write(response.text)
                else:
                    print(f"\n❌ Not found: {uniprot_id} (HTTP {response.status_code})")

            except requests.RequestException as e:
                print(f"\n⚠️ Error retrieving {uniprot_id}: {e}")

            pbar.update(1)  # Update the progress bar
            time.sleep(0.5)  # Prevent overloading UniProt's servers

    print(f"✅ Finished retrieving sequences. Saved to {output_file}")


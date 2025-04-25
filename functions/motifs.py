import matplotlib.pyplot as plt
import biotite.application.clustalo as clustalo
from biotite.application.clustalo import ClustalOmegaApp
import biotite.sequence as seq
import biotite.sequence.graphics as graphics
import numpy as np

from biotite.sequence import ProteinSequence
from biotite.sequence.align import Alignment

from funct.fasta_handling import plot_msa, write_fasta_file, fasta_to_biotite_seq



def motif_logoplot(alignment, size=(4.5, 4.0), axis_pos=None, motif_start=0, png_out=None):
    """
    fig = plt.figure(figsize=size)
    ax = fig.add_subplot(111)
    graphics.plot_alignment_similarity_based(
        ax, alignment[:, :94], labels=sources[:94], symbols_per_line=len(alignment)
    )
    # Source names in italic
    ax.set_yticklabels(ax.get_yticklabels(), fontdict={"fontstyle": "italic"})
    fig.tight_layout()
    """
    profile = seq.SequenceProfile.from_alignment(alignment)

    print("Consensus sequence:")
    consensus = profile.to_consensus()
    print(consensus)

    fig = plt.figure(figsize=size)
    ax = fig.add_subplot(111)
    graphics.plot_sequence_logo(ax, profile, scheme="flower")


    if axis_pos:

        ax.set_xticks(list(range(len(axis_pos))))  # Set correct number of ticks
        ax.set_xticklabels(axis_pos, fontsize=7)  # Assign labels with correct length


    # Increase tick label size
    ax.tick_params(axis='x', labelsize=7)  # Adjust x-axis tick labels
    ax.tick_params(axis='y', labelsize=7)  # Adjust y-axis tick labels
    ax.set_xlabel("Residue Position", fontweight="bold", fontsize=9)
    ax.set_ylabel("Bits", fontweight="bold", fontsize=9)
    # Only show left and bottom spine
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    fig.tight_layout()
    if png_out:
        plt.savefig(png_out, dpi=300)#, dpi=150, bbox_inches='tight')

    #plt.show()
    plt.close()
    return profile.to_consensus()

def get_subsequence(protein_sequence, start, end):
    """
    Extract a subsequence from a ProteinSequence object.

    Parameters:
        protein_sequence (ProteinSequence): The original protein sequence.
        start (int): Start position (0-based index).
        end (int): End position (exclusive, 0-based index).

    Returns:
        ProteinSequence: The subsequence from start to end positions.
    """
    # Validate input
    if not isinstance(protein_sequence, ProteinSequence):
        raise TypeError("Input must be a ProteinSequence object.")
    if start < 0 or end > len(protein_sequence) or start >= end:
        raise ValueError("Invalid start or end positions.")

    # Extract the subsequence
    subsequence = protein_sequence[start:end]
    return ProteinSequence(subsequence)


def extract_subsequences_from_trace(aligned_sequences, alignment_trace, alignment_ids, index_range):
    """
    Extracts subsequences from a list of aligned sequences based on a given index range in the alignment trace.

    Parameters:
        aligned_sequences (list of str or np.ndarray): List of aligned sequences (including gaps).
        alignment_trace (np.ndarray): Alignment trace array (shape: [num_seqs, alignment_length]).
        index_range (tuple): A tuple (start, end) defining the range in the trace.

    Returns:
        list of str: Extracted subsequences for each sequence in the alignment.
    """
    start_pos, end_pos = index_range


    alignment_trace = np.array(alignment_trace)

    start_idx = 0
    end_idx = 0
    for i in range(alignment_trace.shape[0]):
        if alignment_trace[i][0] == start_pos:
            start_idx = i
        if alignment_trace[i][0] == end_pos:
            end_idx = i


    # Validate index range
    if start_idx < 0 or end_idx > alignment_trace.shape[0] or start_idx >= end_idx:
        print(f"Start: {start_idx}; End: {end_idx}; Shape: {alignment_trace.shape[0]}")
        raise ValueError("Invalid index range provided.")

    print(start_idx, end_idx)
    motif_ids = []
    motif_starts =[]
    subsequences = []

    # Extract the relevant columns from the alignment trace
    selected_start_positions = alignment_trace[start_idx, : ] #
    selected_end_positions = alignment_trace[end_idx, :]

    #print(selected_positions)
    if len(selected_start_positions) == len(selected_end_positions):
        for i in range(len(selected_start_positions)):

            if (selected_start_positions[i] != -1):
                seq_start_pos = selected_start_positions[i]
                seq_end_pos = selected_end_positions[i]
                protein_seq = ProteinSequence(aligned_sequences[i])

                length = seq_end_pos - seq_start_pos

                if (length > 0) :
                    subseq = get_subsequence(protein_seq, seq_start_pos, seq_end_pos)
                    if len(subseq) <= length:
                        subsequences.append(subseq)
                        motif_ids.append(alignment_ids[i])
                        motif_starts.append(int(selected_start_positions[i]))

    return subsequences, motif_ids, motif_starts



def extract_subsequences_from_aln(aligned_sequences, alignment_trace, alignment_ids, index_range):
    """
    Extracts subsequences from a list of aligned sequences based on a given index range in the alignment trace.

    Parameters:
        aligned_sequences (list of str or np.ndarray): List of aligned sequences (including gaps).
        alignment_trace (np.ndarray): Alignment trace array (shape: [num_seqs, alignment_length]).
        index_range (tuple): A tuple (start, end) defining the range in the trace.

    Returns:
        list of str: Extracted subsequences for each sequence in the alignment.
    """
    start_pos, end_pos = index_range


    alignment_trace = np.array(alignment_trace)

    start_idx = 0
    end_idx = 0
    for i in range(alignment_trace.shape[0]):
        if i == start_pos:
            start_idx = i
        if i == end_pos:
            end_idx = i


    # Validate index range
    if start_idx < 0 or end_idx > alignment_trace.shape[0] or start_idx >= end_idx:
        print(f"Start: {start_idx}; End: {end_idx}; Shape: {alignment_trace.shape[0]}")
        raise ValueError("Invalid index range provided.")

    print(start_idx, end_idx)
    motif_ids = []
    motif_starts =[]
    subsequences = []

    # Extract the relevant columns from the alignment trace
    selected_start_positions = alignment_trace[start_idx, : ] #
    selected_end_positions = alignment_trace[end_idx, :]

    #print(selected_positions)
    if len(selected_start_positions) == len(selected_end_positions):
        for i in range(len(selected_start_positions)):

            if (selected_start_positions[i] != -1):
                seq_start_pos = selected_start_positions[i]
                seq_end_pos = selected_end_positions[i]
                protein_seq = ProteinSequence(aligned_sequences[i])

                length = seq_end_pos - seq_start_pos

                if (length > 0) :
                    subseq = get_subsequence(protein_seq, seq_start_pos, seq_end_pos)
                    if len(subseq) <= length:
                        subsequences.append(subseq)
                        motif_ids.append(alignment_ids[i])
                        motif_starts.append(int(selected_start_positions[i]))

    return subsequences, motif_ids, motif_starts


def extract_last_subsequences(aligned_sequences, alignment_trace, alignment_ids, last=20):
    """
    Extracts the last `n` residues of each protein sequence based on the alignment trace.

    Parameters:
        aligned_sequences (list of ProteinSequence): List of ProteinSequence objects.
        alignment_trace (np.ndarray): Alignment trace array (shape: [num_seqs, alignment_length]).
        alignment_ids (list): List of sequence identifiers.
        last (int): Number of residues to extract from the end of each sequence (default is 20).

    Returns:
        list of ProteinSequence: Extracted subsequences for each sequence in the alignment.
        list: Corresponding IDs of the sequences.
    """
    subsequences = []
    motif_ids = []
    motif_starts = []

    # Loop through each sequence in the alignment
    for i, protein_seq in enumerate(aligned_sequences):
        # Extract the position in the alignment trace
        end_pos = len(protein_seq)  # Last position of the sequence
        start_pos = max(0, end_pos - last)  # Start position for the last `n` residues

        # Validate that the sequence is a ProteinSequence object
        if not isinstance(protein_seq, ProteinSequence):
            raise TypeError(f"Sequence at index {i} is not of type ProteinSequence.")

        # Get the subsequence using the helper function
        subseq = get_subsequence(protein_seq, start_pos, end_pos)
        subsequences.append(subseq)
        motif_ids.append(alignment_ids[i])
        motif_starts.append(start_pos)

    return subsequences, motif_ids, motif_starts

def msa_of_motif(file):
    """
    Performs multiple sequence alignment (MSA) for a given motif using Clustal Omega.

    Parameters:
        name (str): Name of the motif.
        file (str): Path to the FASTA file containing sequences.
        starts (dict): Dictionary containing start positions of sequences.

    Returns:
        tuple: (MSA alignment object, dictionary of sequences, list of sequences,
                list of sequence IDs, list of start positions).
    """
    hit_dict, hit_seqs, hit_ids, hit_starts = fasta_to_biotite_seq(file)#, starts)

    aln = clustalo.ClustalOmegaApp.align(hit_seqs, bin_path="C:/Users/c-fih/anaconda3/envs/clustalo/clustalo-122-win64/clustal-omega-1.2.2-win64/clustalo.exe")

    return aln, hit_dict, hit_seqs, hit_ids, hit_starts




def msa_of_motif_strict(file):
    """
    Performs multiple sequence alignment (MSA) for a given motif using Clustal Omega.

    Parameters:
        name (str): Name of the motif.
        file (str): Path to the FASTA file containing sequences.
        starts (dict): Dictionary containing start positions of sequences.

    Returns:
        tuple: (MSA alignment object, dictionary of sequences, list of sequences,
                list of sequence IDs, list of start positions).
    """
    hit_dict, hit_seqs, hit_ids, hit_starts = fasta_to_biotite_seq(file)#, starts)

    app = ClustalOmegaApp(hit_seqs, bin_path="C:/Users/c-fih/anaconda3/envs/clustalo/clustalo-122-win64/clustal-omega-1.2.2-win64/clustalo.exe")
    #app.add_additional_options("")
    app.start()
    app.join()
    aln = app.get_alignment()
    return aln, hit_dict, hit_seqs, hit_ids, hit_starts

def motif_analysis_from_aln(aln, aln_ids, motif_name, motif_file, motif_pos, readonly=True):
    """
    Extracts a motif from an alignment and performs multiple sequence alignment (MSA) analysis.

    Parameters:
        aln (Alignment): Input sequence alignment object.
        aln_ids (list): List of sequence IDs in the alignment.
        motif_name (str): Name of the motif.
        motif_file (str): Path to save the motif sequences in FASTA format.
        motif_pos (tuple): Start and end positions of the motif within the alignment.
        readonly (bool, optional): If False, only reads motifs from existing FASTA file.
                Default is True, extracts motifs from alignment and saves them into FASTA file.

    Returns:
        tuple: (MSA alignment of the motif, dictionary of motif sequences, list of motif sequences,
                list of motif sequence IDs, list of motif start positions).
    """

    print(f"---------------{motif_name} motif:")
    if not readonly:
        motif, motif_ids, motif_starts = extract_subsequences_from_trace(aln.sequences, aln.trace, aln_ids, motif_pos)
        write_fasta_file(motif, motif_ids, f"{motif_file}.fasta", description=motif_name)

    motif_aln, motif_dict, motif_seqs, motif_ids, motif_starts = msa_of_motif(f"{motif_file}.fasta")#, motif_starts)
    plot_msa(motif_aln, motif_ids, motif_starts, title=f"MSA of {motif_name} motif", png_out=f"{motif_file}.png")

    return motif_aln, motif_dict, motif_seqs, motif_ids, motif_starts



def c_terminus_analysis_from_aln(aln, aln_ids, motif_name, motif_file, cterm_len=20, readonly=True):
    """
    Extracts and analyzes the C-terminal region of aligned sequences.

    Parameters:
        aln (Alignment): Input sequence alignment object.
        aln_ids (list): List of sequence IDs in the alignment.
        motif_name (str): Name of the motif.
        motif_file (str): Path to save extracted sequences in FASTA format.
        cterm_len (int, optional): Length of the C-terminal region to extract. Default is 20.
        readonly (bool, optional): If False, extracts and saves the C-terminal regions in a FASTA file. Default is True.
        readonly (bool, optional): If False, only reads C-terminal regions from existing FASTA file.
                Default is True, extracts motifs from alignment and saves them into FASTA file.

    Returns:
        tuple: (MSA alignment of C-terminal region, dictionary of extracted sequences,
                list of extracted sequences, list of sequence IDs, list of sequence start positions).
    """
    print(f"---------------{motif_name} motif:")
    if not readonly:
        #motif, motif_ids = extract_c_terminus_from_trace(aln.sequences, aln.trace, aln_ids, last=cterm_len)
        motif, motif_ids, motif_starts = extract_last_subsequences(aln.sequences, aln.trace, aln_ids, last=cterm_len)
        #for seq, motif_id in zip(motif, motif_ids):
            #print(f">{motif_id}\n{seq}")
        write_fasta_file(motif, motif_ids, f"{motif_file}.fasta", description=motif_name)

    motif_aln, motif_dict, motif_seqs, motif_ids, motif_starts = msa_of_motif(f"{motif_file}.fasta")#, f"{motif_file}.fasta", motif_starts)
    plot_msa(motif_aln, motif_ids, motif_starts, title=f"MSA of {motif_name} motif", png_out=f"{motif_file}.png")

    return motif_aln, motif_dict, motif_seqs, motif_ids, motif_starts

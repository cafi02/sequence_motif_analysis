import biotite

import biotite.sequence.graphics as graphics
import matplotlib.pyplot as plt


from Bio import SeqIO
from Bio import AlignIO
from biotite.sequence import ProteinSequence
from biotite.sequence.align import Alignment
import numpy as np



def add_fasta2fasta(source1_fasta, source2_fasta, target_fasta):
    """
    Appends the contents of source_fasta to the end of target_fasta.

    Parameters:
    source_fasta (str): Path to the source FASTA file to append.
    target_fasta (str): Path to the target FASTA file to which sequences will be added.
    """
    try:
        with open(source1_fasta, "r") as src1, open(source2_fasta, "r") as src2, open(target_fasta, "w") as tgt:
            tgt.write(src1.read())
            tgt.write("\n")
            tgt.write(src2.read())
        print(f"✅ Successfully appended {source1_fasta} and {source2_fasta} to:\n{target_fasta}")
    except Exception as e:
        print(f"Error: {e}")



def fasta_to_biotite_seq(file_path, full_headers=False, starts=None):
    """
    Converts sequences from a FASTA file to Biotite ProteinSequence objects.

    Args:
    file_path (str): The path to the FASTA file.

    Returns:
    tuple: A dictionary of sequences from the FASTA file, a list of Biotite ProteinSequence objects,
           a list of sequence IDs, and a list of sequence start positions.
    """

    blast_dict = SeqIO.to_dict(SeqIO.parse(file_path, "fasta"))
    hit_seqs = []
    hit_ids = []
    hit_starts = []

    s = 0
    for hit in blast_dict.values():
        seq = str(hit.seq).strip("-").replace("-", "")

        if ("U" not in seq) :
            sequence = biotite.sequence.ProteinSequence(seq)  # convert into biotite ProteinSequence object
            hit_seqs.append(sequence)
            if full_headers:
                try:
                    with open(file_path, 'r') as file:
                        for line in file:
                            if hit.id in line:
                                id = line.strip('>tr|').split(' OX')[0]
                                hit_ids.append(id)
                except FileNotFoundError:
                    id = hit.id
                    hit_ids.append(id)

            else:
                id = hit.id
                hit_ids.append(id)

            if starts:
                hit_starts.append(starts[s])
                s += 1
            else:
                hit_starts.append(1)


    print(f"Number of Hits in {file_path}:\t{len(hit_ids)}")

    return blast_dict, hit_seqs, hit_ids, hit_starts



def write_fasta_file(protein_sequences, headers, output_file, description="Sequence motif"):
    """
    Write a list of ProteinSequence objects and corresponding headers to a FASTA file.

    Parameters:
        protein_sequences (list): List of ProteinSequence objects.
        headers (list): List of headers (without '>').
        output_file (str): Path to the output FASTA file.

    Raises:
        ValueError: If the lengths of protein_sequences and headers do not match.
    """
    if len(protein_sequences) != len(headers):
        raise ValueError("The number of sequences must match the number of headers.")

    if not all(isinstance(seq, ProteinSequence) for seq in protein_sequences):
        raise TypeError("All items in protein_sequences must be of type ProteinSequence.")

    with open(output_file, 'w') as fasta_file:
        for header, protein_seq in zip(headers, protein_sequences):
            fasta_file.write(f">{header}\n{str(protein_seq)}\n")

    print(f"FASTA file for {description} saved as: {output_file}")



def update_fasta_headers(motif_fasta, full_fasta):
    # Output file name
    output_file = motif_fasta.split(".")[0] + "_descr.fasta"

    # Parse the FASTA files
    motif_records = SeqIO.parse(motif_fasta, "fasta")
    full_records = SeqIO.to_dict(SeqIO.parse(full_fasta, "fasta"))

    with open(output_file, "w") as output_handle:
        for motif_record in motif_records:
            motif_id = motif_record.id
            if motif_id in full_records:
                full_header = full_records[motif_id].description
                # Write the full header from full_fasta
                output_handle.write(f">{full_header}\n")
                # Write the sequence from motif_fasta
                output_handle.write(str(motif_record.seq) + "\n")

    print(f"Output saved to {output_file}")


def plot_msa(alignment, hit_ids, hit_starts, size=(20, 200), title="MSA", png_out=None):
    # Adapted from code source: Patrick Kunzmann
    # License: BSD 3 clause
    """
    Generates and plots a multiple sequence alignment (MSA) visualization.
    If subset=True, only the first 10 and last 10 sequences are shown in separate subplots.

    Parameters:
        alignment (Alignment): Multiple sequence alignment object.
        hit_ids (list): List of sequence IDs.
        hit_starts (list): List of start positions for sequences.
        size (tuple, optional): Figure size. Default is (20, 200).
        title (str, optional): Title for the plot. Default is "MSA".
        png_out (str, optional): File path to save the plot. Default is None.
        subset (bool, optional): If True, only plots the first 10 and last 10 sequences in separate subplots.

    Returns:
        None: Displays the plot or saves it if png_out is provided.
    """
    number_functions = []
    trace_length = len(alignment.trace[:])
    symbols_per_line = 100 if trace_length > 100 else trace_length

    for start in hit_starts:
        def some_func(x, start=start):
            return x + start

        number_functions.append(some_func)

    ids = [hit.split("|")[1] for hit in hit_ids]

    # Full alignment plot
    fig, ax = plt.subplots(figsize=size)
    graphics.plot_alignment_type_based(
        ax,
        alignment,
        symbols_per_line=symbols_per_line,
        labels=ids,
        symbol_size=7,
        number_size=7,
        label_size=7,
        show_numbers=True,
        number_functions=number_functions,
        color_scheme="flower",
    )
    ax.set_title(title)

    # Save the plot if output path is provided
    if png_out:
        plt.savefig(png_out)#, dpi=1200)
    plt.close(fig)


def read_alignment_from_fasta(fasta_file):
    """
    Reads an aligned FASTA file and converts it into a Biotite Alignment object.

    Parameters:
    fasta_file (str): Path to the aligned FASTA file.

    Returns:
    biotite.sequence.align.Alignment: Alignment object containing sequences with gaps.
    """
    align = AlignIO.read(fasta_file, "fasta")
    trace = []
    sequences = []
    for record in align:
        seq_trace = []
        sequence = record.seq
        sequences.append(ProteinSequence(sequence.replace("-", "")))

        pos = 0
        for i in range(len(sequence)):
            if not sequence[i] == "-":
                seq_trace.append(pos)
                pos += 1
            else:
                seq_trace.append(-1)

        trace.append(seq_trace)

    trace = np.transpose(np.array(trace))

    alignment = Alignment(sequences, trace)

    return alignment
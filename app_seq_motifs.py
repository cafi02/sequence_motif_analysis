from funct.fasta_handling import fasta_to_biotite_seq, plot_msa, read_alignment_from_fasta
from funct.coevolution import mutual_information_zscore, plot_coev
from funct.motifs import motif_logoplot
import numpy as np


# FUNCTION
def motif_analysis(motif_name, file, motif_out_file, positions=None):
    # Read MSA from FASTA file (e.g. extracted from Jalviewer
    motif_aln = read_alignment_from_fasta(file)
    motif_dict, motif_seqs, motif_ids, motif_starts = fasta_to_biotite_seq(file)

    # Plot MSA
    plot_msa(motif_aln, motif_ids, motif_starts, title=f"MSA of {motif_name}", png_out=f"{motif_out_file}_msa.png")
    print(f"✅ Finished MSA for: {motif_name}")

    # Logoplot
    #consensus = motif_logoplot(motif_aln, size=(8, 1.7), axis_pos=positions, png_out=f"{motif_out_file}_logoplot.png")
    consensus = motif_logoplot(motif_aln, size=(3.5, 1.7), axis_pos=positions, png_out=f"{motif_out_file}_logoplot.png")
    print(f"✅ Finished Logoplot for: {motif_name}")

    # Coevolution
    motif_aln = motif_aln[motif_aln.trace[:, 0] != -1]  # Remove alignment columns that have a gap in the sequence
    motif_mi = mutual_information_zscore(motif_aln)
    np.savetxt(f"{motif_out_file}_coev.txt", motif_mi, delimiter="\t", fmt="%.6f")
    plot_coev(motif_mi,  title=f"{motif_name} - Coevolution", axis_seq=consensus, axis_pos=positions,
              png_out=f"{motif_out_file}_coev.png")
    print(f"✅ Finished Coevolution for: {motif_name}")



# EDIT INPUT PARAMS
motif_name = "Water Channel"
motif_infile = "data/clustalo/phmmer_rp15_key1_E-50_jalview_motif1.fasta"
motif_outfile = f"data/seq_motifs/jalview/rp15_key1_E-50_jalview_motif1"
motif_pos = [0, 8, 9, 10, 11, 12, 13, 14, 25, 26, 27, 28] # no continious positions

# APPLICATION
motif_analysis(motif_name, motif_infile, motif_outfile, positions=motif_pos)

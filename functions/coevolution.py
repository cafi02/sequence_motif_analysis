# Adapted from code source: Patrick Kunzmann
# License: BSD 3 clause

import warnings
import biotite.sequence.align as align
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import biotite
from mpl_toolkits.axes_grid1 import make_axes_locatable



def mutual_information_zscore(alignment, n_shuffle=100):
    """
    Code source: Patrick Kunzmann (License: BSD 3 clause)
    Computes the mutual information Z-score matrix for a sequence alignment.

    Parameters:
        alignment (Alignment): Multiple sequence alignment object.
        n_shuffle (int, optional): Number of shuffles for computing random MI distributions. Default is 100.

    Returns:
        numpy.ndarray: Z-score matrix representing coevolutionary relationships between residues.
    """
    codes = align.get_codes(alignment).T
    alph = alignment.sequences[0].alphabet

    mi = _mutual_information(alignment, codes, alph)
    np.random.seed(0)
    random_mi = [None] * n_shuffle
    for i in range(n_shuffle):
        shuffled_codes = _shuffle(codes)
        random_mi[i] = _mutual_information(alignment, shuffled_codes, alph)
    random_mi = np.stack(random_mi)
    mean = np.mean(random_mi, axis=0)
    std = np.std(random_mi, axis=0)
    z_score = (mi - mean) / std
    return z_score


def _shuffle(codes):
    """
    Code source: Patrick Kunzmann (License: BSD 3 clause)
    Randomly shuffles each column of an alignment matrix.

    Parameters:
        codes (numpy.ndarray): Matrix representation of sequence alignment.

    Returns:
        numpy.ndarray: Shuffled matrix.
    """

    shuffled_codes = codes.copy()
    # Shuffle each alignment column
    for i in range(len(shuffled_codes)):
        np.random.shuffle(shuffled_codes[i])
    return shuffled_codes


def _mutual_information(alignment, codes, alph):
    """
    Code source: Patrick Kunzmann (License: BSD 3 clause)
    Computes the mutual information matrix for a given sequence alignment.

    Parameters:
        alignment (Alignment): Multiple sequence alignment object.
        codes (numpy.ndarray): Matrix representation of sequence alignment.
        alph (Alphabet): Alphabet object representing sequence characters.

    Returns:
        numpy.ndarray: Mutual information matrix capturing coevolutionary signals.
    """

    mi = np.zeros((len(alignment), len(alignment)))
    # Iterate over all columns to choose first column
    for i in range(codes.shape[0]):
        # Iterate over all columns to choose second column
        for j in range(codes.shape[0]):
            nrows = 0
            marginal_counts_i = np.zeros(len(alph), dtype=int)
            marginal_counts_j = np.zeros(len(alph), dtype=int)
            combined_counts = np.zeros((len(alph), len(alph)), dtype=int)
            # Iterate over all symbols in both columns
            for k in range(codes.shape[1]):
                # Skip rows where either column has a gap
                if codes[i, k] != -1 and codes[j, k] != -1:
                    marginal_counts_i[codes[i, k]] += 1
                    marginal_counts_j[codes[j, k]] += 1
                    combined_counts[codes[i, k], codes[j, k]] += 1
                    nrows += 1
            marginal_probs_i = marginal_counts_i / nrows
            marginal_probs_j = marginal_counts_j / nrows
            combined_probs = combined_counts / nrows

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                mi_before_sum = (
                    combined_probs
                    * np.log2(
                        combined_probs
                        / (
                            marginal_probs_i[:, np.newaxis]
                            * marginal_probs_j[np.newaxis, :]
                        )
                    )
                ).flatten()
            mi[i, j] = np.sum(mi_before_sum[~np.isnan(mi_before_sum)])
    return mi


def combine_string_and_list(int_list, char_string):
    if len(int_list) != len(char_string):
        raise ValueError("Input positions and consensus are not of the same length!")

    return [f"{char}{num}" for char, num in zip(char_string, int_list)]


def plot_coev(mi, size=None, motif_start=0, axis_seq=None, axis_pos=None, title="Residue Coevolution", png_out=None):
    """
    Adapted from code source: Patrick Kunzmann (License: BSD 3 clause)
    Plots a coevolution matrix based on mutual information Z-scores.

    Parameters:
        mi (numpy.ndarray): MI-Z-score matrix.
        title (str, optional): Title of the plot. Default is "Residue Coevolution".
        png_out (str, optional): File path to save the plot. Default is None.

    Returns:
        None: Displays or saves the coevolution matrix plot.
    """
    if axis_seq and not size:
        i = (len(axis_seq) * 0.15) + 1.2
        size = (i, i)
    elif not size:
        size = (5, 5)

    # Create the color map for the plot
    color = colors.to_rgb(biotite.colors["dimorange"])
    cmap_val = np.stack(
        [
            np.interp(np.linspace(0, 1, 100), [0, 1], [1, color[i]])
            for i in range(len(color))
        ]
    ).transpose()
    cmap = colors.ListedColormap(cmap_val)

    # Create figure and axis
    fig, ax = plt.subplots(figsize=size)

    # Create a divider to ensure the colorbar has the same height as the plot
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)  # Adjust size and padding

    # Plot the MI-Z-Score matrix
    im = ax.pcolormesh(mi, cmap=cmap)
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("MI-Z-score", fontsize=11)

    ax.set_aspect("equal")
    if axis_seq and axis_pos:
        ax.set_xlabel("Consensus sequence", fontsize=9, fontweight="bold")
        ax.set_ylabel("Consensus sequence", fontsize=9, fontweight="bold")
    else:
        ax.set_xlabel("Residue position", fontsize=9, fontweight="bold")
        ax.set_ylabel("Residue position", fontsize=9, fontweight="bold")
    ax.set_title(title, fontsize=11, pad=18)

    if axis_seq and axis_pos and len(axis_seq)==mi.shape[0]:
        # Set x and y ticks at the center of each corresponding heat map cell
        tick_positions = np.arange(0.5, mi.shape[0])  # Centered ticks
        #tick_labels = [aa for aa in axis_seq]  # consensus sequence
        print(axis_seq, axis_pos)
        tick_labels = combine_string_and_list(axis_pos[1:], axis_seq)
        print(tick_labels)
    else:
        # Set x and y ticks at the center of each corresponding heat map cell
        steps = int(round(mi.shape[0] / 5, 0))
        tick_positions = np.arange(0.5, mi.shape[0], steps)  # Centered ticks
        tick_labels = np.arange(motif_start, motif_start + mi.shape[0], steps) + 1  # Increase labels by 1

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontsize=7, rotation=90)  # Ensure integer labels
    ax.set_yticks(tick_positions)
    ax.set_yticklabels(tick_labels, fontsize=7)  # Ensure integer labels

    fig.tight_layout()

    if png_out:
        plt.savefig(png_out, dpi=500)

    #plt.show()
    plt.close()



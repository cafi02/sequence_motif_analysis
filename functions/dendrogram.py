# Adapted from code source: Patrick Kunzmann
# License: BSD 3 clause

import matplotlib.pyplot as plt

import biotite.sequence.align as align
import biotite.sequence.phylo as phylo
import biotite.sequence.graphics as graphics




def abbreviate(species):
    """
    Abbreviates a species name by using the first letter of the genus and the full species name.

    Parameters:
        species (str): Full species name.

    Returns:
        str: Abbreviated species name in the format "G. species", or the original name if it cannot be abbreviated.
    """
    # Remove possible brackets
    species = species.replace("[", "").replace("]", "").strip()
    splitted_species = species.split()

    # Ensure there are at least two parts
    if len(splitted_species) < 2:
        return species  # Return the original name if abbreviation is not possible

    return f"{splitted_species[0][0]}. {splitted_species[1]}"


def get_hit_description(list_of_ids):
    short_ids = []

    for i in list_of_ids:
        hit = i.split('|')[0]
        enzyme = ((i.split("OS=")[0]).split("|")[-1]).strip().split(" ", 1)[
            -1]  # f'{i.split(" ")[1]} {i.split(" ")[2]}'
        organism = f'{(i.split("OS=")[-1]).split(" ")[0]} {(i.split("OS=")[-1]).split(" ")[1]}'
        organism = abbreviate(organism)
        description = f'{hit} | {enzyme} | {organism}'

        short_ids.append(description)

    return short_ids

def dendrogram_from_aln(alignment, ids, title="Sequence deviation", png_out=None):
    # Get distance matrix:
    distances = 1 - align.get_pairwise_sequence_identity(alignment, mode="shortest")
    tree = phylo.upgma(distances)

    ### plot the tree
    fig, ax = plt.subplots(1, 1, figsize=(8, 15))
    graphics.plot_dendrogram(
        ax,
        tree,
        orientation="top",
        labels=list(ids),
        show_distance=False,
        linewidth=2,
    )
    ax.grid(False)
    ax.set_yticks([])

    for label in ax.get_xticklabels():
        label.set_rotation(-90)
        label.set_horizontalalignment("center")

    plt.subplots_adjust(left=0.01, right=0.99, bottom=0.5, top=0.9)

    # distance indicator
    indicator_len = 0.1
    indicator_start = (0,0)

    indicator_stop = (indicator_start[0], indicator_start[1] + indicator_len)
    indicator_center = (
        (indicator_start[0] + 1,
        ((indicator_start[1] + indicator_stop[1]) / 2))
    )
    ax.annotate(
        "",
        xy=indicator_start,
        xytext=indicator_stop,
        xycoords="data",
        textcoords="data",
        arrowprops={"arrowstyle": "|-|", "linewidth": 2},
    )
    ax.annotate(
        f"{int(indicator_len * 100)} %", xy=indicator_center, ha="center", va="center"
    )
    ax.set_title(title)
    if png_out:
        plt.savefig(png_out, dpi=300, bbox_inches="tight")

    #plt.show()
    plt.close()

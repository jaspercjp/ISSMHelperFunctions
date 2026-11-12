import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import matplotlib.colors as mcolors

def plot_field(x, y, v, figsize=(10, 5), ax=None, cmap=plt.get_cmap("viridis"), clabel='v', vmin=-100, vmax=100, block_out=False):
    """
    Plot scalar field v at coordinates (x, y) using triangulation.
    """
    triang = tri.Triangulation(x, y)
    if ax is None:
        plt.figure(figsize=figsize)

    if ax is None:
        tpc = plt.tripcolor(triang, v, shading="flat", cmap=cmap, vmin=vmin, vmax=vmax)
    else:
        if block_out:
            # Create a white colormap
            white_cmap = mcolors.ListedColormap(['red'])

            # Plot with white colormap
            tpc = ax.tripcolor(triang, v, shading="flat", cmap=white_cmap, vmin=vmin, vmax=vmax)
        else:
            tpc = ax.tripcolor(triang, v, shading="flat", cmap=cmap, vmin=vmin, vmax=vmax)

    return tpc
    # plt.colorbar(tpc, label=clabel)
    # plt.xlabel("x")
    # plt.ylabel("y")
    # plt.axis("equal")
    # plt.xlim(min(x), max(x))
    # plt.ylim(min(y), max(y))

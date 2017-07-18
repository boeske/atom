"""
ATOM: Plot function
Version 2.0
"""

# Numpy
import numpy as np
from numpy import unravel_index

# Scipy
import scipy as sp
from scipy import ndimage, interpolate
from scipy.interpolate import interp2d

# Matplotlib
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib import rc

from source.messages import *
from source.conf import *


def atom_plot_new(topo, stage, number):
    ''' The new default topography plotting function '''

    # Print a welcome message
    message_warning("ATOM plot: stage " + stage + ", (plot #" + str(number) + ")")

    # Initialize the figure
    fig = plt.figure()

    # Set title of plot
    plt.title(case_name + ": " + stage)

    # pick the desired colormap, sensible levels, and define a normalization
    # instance which takes data values and translates those into levels.
    levels = MaxNLocator(nbins=100).tick_values(topo[stage]["data"].min(),
                                                topo[stage]["data"].max())
    cmap = plt.get_cmap('viridis')
    # cmap = plt.get_cmap('terrain')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    y, x = np.mgrid[slice(0, topo[stage]["Ly"] - 1, topo[stage]["dx"]),
                    slice(0, topo[stage]["Lx"] - 1, topo[stage]["dx"])]

    # Make the plot
    plt.pcolormesh(x, y, topo[stage]["data"], cmap=cmap, norm=norm)

    # Setup x/y axis
    plt.axes().set_aspect('equal')
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    plt.axes().set_xlabel('x (m)')
    plt.axes().set_ylabel('y (m)')

    # Draw a colorbar
    cbar = plt.colorbar(orientation="horizontal")
    cbar.set_label('Height (m)')

    # Mark wind turbines with green stars (UTM values)
    for i in range(1, n_turbines + 1):
        ID = str(i)
        plt.plot(turbines[stage][ID]["utm_x"] - topo["raw"]["origin_x"],
                 turbines[stage][ID]["utm_y"] - topo["raw"]["origin_y"], 'g2')

    # Mark wind turbines with red stars (After rotation this is needed)
    for i in range(1, n_turbines + 1):
        ID = str(i)
        plt.plot(turbines[stage][ID]["i"] * topo[stage]["dx"],
                 turbines[stage][ID]["j"] * topo[stage]["dx"], 'r2')

    # Mark wind turbines with yellow stars
    for i in range(1, n_masts + 1):
        ID = str(i)
        plt.plot(masts[stage][ID]["utm_x"] - topo["raw"]["origin_x"],
                 masts[stage][ID]["utm_y"] - topo["raw"]["origin_y"], 'gx')

    # Mark masts with blue stars (After rotation this is needed)
    for i in range(1, n_masts + 1):
        ID = str(i)
        plt.plot(masts[stage][ID]["i"] * topo[stage]["dx"],
                 masts[stage][ID]["j"] * topo[stage]["dx"], 'bx')

        # Temporary: Plot lower left coordinate of final PALM domain
        plt.plot(crosssections[stage][ID]["origin_utm_x"] -
                 topo["raw"]["origin_x"],
                 crosssections[stage][ID]["origin_utm_y"] -
                 topo["raw"]["origin_y"], 'bs')

        if nesting:
            if stage == "PALM":
                plt.plot(crosssections[stage]["1"]["nest_origin_i"] * topo[stage]["dx"],
                         crosssections[stage]["1"]["nest_origin_j"] * topo[stage]["dx"], 'rd')

    # Save the figure
    fig.savefig(output_dir + str(number) + "_" + stage + ".png", dpi=dpi_value)

    # Delete the figure to free memory
    plt.close('all')

#!/usr/bin/python3
"""
ATOM - the Almighty TOpography Modifier
Version 2.0
"""


# ------------------------------------------------------------------------------
#                          Import required libraries
# ------------------------------------------------------------------------------

# ATOM variables and functions
from source.conf import *
from source.messages import *
from source.functions import *
from source.plot import *

# Numpy
import numpy as np
from numpy import unravel_index

# Scipy
import scipy as sp
from scipy import ndimage, interpolate
from scipy.interpolate import interp2d

# Others
import os
import os.path
import math
from io import StringIO
from os import mkdir
import shutil
from shutil import copyfile
import sys
import pprint


# ------------------------------------------------------------------------------
#                               Start ATOM
# ------------------------------------------------------------------------------
startmessage()

# Create output directory
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Copy atom.ini to output-directory
ini_path = output_dir + "atom.ini"
copyfile("atom.ini", ini_path)

# Create PALM-INPUT directory
palm_input_dir = "output/" + case_name + "/INPUT"
if not os.path.exists(palm_input_dir):
    os.makedirs(palm_input_dir)


# ------------------------------------------------------------------------------
#                      Print out main settings
# ------------------------------------------------------------------------------
message("Settings")
print("  Input file: " + str(inputfile))
print("  p3d-template: " + str(p3d_template))
print("  Grid spacing: " + str(dx) + " m")
print("  Wind direction: " + str(wdir) + "°")
print("  Rotation angle: " + str(angle) + "°")
print("  Number of turbines: " + str(n_turbines))
print("  Number of masts: " + str(n_masts))

if nesting:
    print("  --------------------")
    print("  Nesting: " + boldgreen("ON"))
    print("  dx_nest: " + str(dx_nest) + " m")
    print("  --------------------")

print("  Output to: output/" + bold(case_name))


# ------------------------------------------------------------------------------
#                             Read ASCII file
# ------------------------------------------------------------------------------
message("Read ASCII file...")

# Read ASCII file
message_info("Read data")
f = open(inputfile)
data = f.read()
f.close()
message_done()

# Generate a numpy data array
message_info("Generate data array")
data_str = StringIO(data)
topo["raw"]["data"] = np.loadtxt(data_str)
message_done()

# Flip upside down (array is initially upside down because of ASCII format)
message_info("Flipping data upside down")
topo["raw"]["data"] = np.flipud(topo["raw"]["data"])
message_done()

# Store basic parameters of raw input topography
topo["raw"]["dx"] = dx_raw
topo["raw"]["nx"] = topo["raw"]["data"].shape[1]
topo["raw"]["ny"] = topo["raw"]["data"].shape[0]
topo["raw"]["Lx"] = topo["raw"]["nx"] * topo["raw"]["dx"]
topo["raw"]["Ly"] = topo["raw"]["ny"] * topo["raw"]["dx"]
topo["raw"]["origin_x"] = origin_x
topo["raw"]["origin_y"] = origin_y
topo["raw"]["center"] = np.array([topo["raw"]["nx"] / 2,
                                  topo["raw"]["ny"] / 2])

# Get the initial turbine coordinates in UTM format
for i in range(1, n_turbines + 1):
    ID = str(i)
    turbines["raw"][ID]["utm_x"] = wka_ij[i - 1, 0]
    turbines["raw"][ID]["utm_y"] = wka_ij[i - 1, 1]
    turbines["raw"][ID]["i"] = (turbines["raw"][ID]["utm_x"] -
                                topo["raw"]["origin_x"]) / topo["raw"]["dx"]
    turbines["raw"][ID]["j"] = (turbines["raw"][ID]["utm_y"] -
                                topo["raw"]["origin_y"]) / topo["raw"]["dx"]

# Get the initial mast coordinates in UTM format
for i in range(1, n_masts + 1):
    ID = str(i)
    masts["raw"][ID]["utm_x"] = masts_ij[i - 1, 0]
    masts["raw"][ID]["utm_y"] = masts_ij[i - 1, 1]
    masts["raw"][ID]["i"] = (masts["raw"][ID]["utm_x"] -
                             topo["raw"]["origin_x"]) / topo["raw"]["dx"]
    masts["raw"][ID]["j"] = (masts["raw"][ID]["utm_y"] -
                             topo["raw"]["origin_y"]) / topo["raw"]["dx"]

# Plot the raw topography
if plotting:
    atom_plot_new(topo, "raw", plot_number)
    plot_number += 1


# ------------------------------------------------------------------------------
#                            Interpolation
# ------------------------------------------------------------------------------
message("Interpolation, dx_raw -> dx")

# Basic parameters of interpolated topography
topo["interpolated"]["dx"] = dx
topo["interpolated"]["nx"] = int(topo["raw"]["Lx"] /
                                 topo["interpolated"]["dx"])
topo["interpolated"]["ny"] = int(topo["raw"]["Ly"] /
                                 topo["interpolated"]["dx"])
topo["interpolated"]["Lx"] = (topo["interpolated"]["nx"] *
                              topo["interpolated"]["dx"])
topo["interpolated"]["Ly"] = (topo["interpolated"]["ny"] *
                              topo["interpolated"]["dx"])
topo["interpolated"]["origin_x"] = topo["raw"]["origin_x"]
topo["interpolated"]["origin_y"] = topo["raw"]["origin_y"]

# Generate grids for interpolation
x_in = np.linspace(0, topo["raw"]["Lx"], topo["raw"]["nx"])
y_in = np.linspace(0, topo["raw"]["Ly"], topo["raw"]["ny"])
x_int = np.linspace(0, topo["interpolated"]["Lx"], topo["interpolated"]["nx"])
y_int = np.linspace(0, topo["interpolated"]["Ly"], topo["interpolated"]["ny"])

message_info('Interpolate...')

# Interpolation function
# It works, but i don't get the syntax here...
int2d = interp2d(x_in, y_in, topo["raw"]["data"], kind='linear')
topo["interpolated"]["data"] = int2d(x_int, y_int)

# Coordinates for turbines and masts do not change from previous stage,
# coordinates in gridpoints however do change
for i in range(1, n_turbines + 1):
    ID = str(i)
    turbines["interpolated"][ID]["utm_x"] = turbines["raw"][ID]["utm_x"]
    turbines["interpolated"][ID]["utm_y"] = turbines["raw"][ID]["utm_y"]
    turbines["interpolated"][ID]["i"] = \
        (turbines["interpolated"][ID]["utm_x"] -
         topo["interpolated"]["origin_x"]) / topo["interpolated"]["dx"]
    turbines["interpolated"][ID]["j"] = \
        (turbines["interpolated"][ID]["utm_y"] -
         topo["interpolated"]["origin_y"]) / topo["interpolated"]["dx"]

for i in range(1, n_masts + 1):
    ID = str(i)
    masts["interpolated"][ID]["utm_x"] = masts["raw"][ID]["utm_x"]
    masts["interpolated"][ID]["utm_y"] = masts["raw"][ID]["utm_y"]
    masts["interpolated"][ID]["i"] = ((masts["interpolated"][ID]["utm_x"] -
                                       topo["interpolated"]["origin_x"]) /
                                      topo["interpolated"]["dx"])
    masts["interpolated"][ID]["j"] = ((masts["interpolated"][ID]["utm_y"] -
                                       topo["interpolated"]["origin_y"]) /
                                      topo["interpolated"]["dx"])

message_done()

# Plot the interpolated topography
if plotting:
    atom_plot_new(topo, "interpolated", plot_number)
    plot_number += 1

# ------------------------------------------------------------------------------
#                   Interpolation for the nested domain
# ------------------------------------------------------------------------------
if nesting:
    message("Interpolation (Nest), dx_raw -> dx_nest")

    # Basic parameters of interpolated topography
    topo["interpolated_nest"]["dx"] = dx_nest
    topo["interpolated_nest"]["nx"] = int(topo["raw"]["Lx"] /
                                          topo["interpolated_nest"]["dx"])
    topo["interpolated_nest"]["ny"] = int(topo["raw"]["Ly"] /
                                          topo["interpolated_nest"]["dx"])
    topo["interpolated_nest"]["Lx"] = (topo["interpolated_nest"]["nx"] *
                                       topo["interpolated_nest"]["dx"])
    topo["interpolated_nest"]["Ly"] = (topo["interpolated_nest"]["ny"] *
                                       topo["interpolated_nest"]["dx"])
    topo["interpolated_nest"]["origin_x"] = topo["raw"]["origin_x"]
    topo["interpolated_nest"]["origin_y"] = topo["raw"]["origin_y"]

    # Generate grids for interpolation
    x_in = np.linspace(0, topo["raw"]["Lx"], topo["raw"]["nx"])
    y_in = np.linspace(0, topo["raw"]["Ly"], topo["raw"]["ny"])
    x_int = np.linspace(0, topo["interpolated_nest"]["Lx"], topo["interpolated_nest"]["nx"])
    y_int = np.linspace(0, topo["interpolated_nest"]["Ly"], topo["interpolated_nest"]["ny"])

    message_info('Interpolate (nest)...')

    # Interpolation function
    # It works, but i don't get the syntax here...
    int2d = interp2d(x_in, y_in, topo["raw"]["data"], kind='linear')
    topo["interpolated_nest"]["data"] = int2d(x_int, y_int)

    # Coordinates for turbines and masts do not change from previous stage,
    # coordinates in gridpoints however do change
    for i in range(1, n_turbines + 1):
        ID = str(i)
        turbines["interpolated_nest"][ID]["utm_x"] = turbines["raw"][ID]["utm_x"]
        turbines["interpolated_nest"][ID]["utm_y"] = turbines["raw"][ID]["utm_y"]
        turbines["interpolated_nest"][ID]["i"] = \
            (turbines["interpolated_nest"][ID]["utm_x"] -
             topo["interpolated_nest"]["origin_x"]) / topo["interpolated_nest"]["dx"]
        turbines["interpolated_nest"][ID]["j"] = \
            (turbines["interpolated_nest"][ID]["utm_y"] -
             topo["interpolated_nest"]["origin_y"]) / topo["interpolated_nest"]["dx"]

    for i in range(1, n_masts + 1):
        ID = str(i)
        masts["interpolated_nest"][ID]["utm_x"] = masts["raw"][ID]["utm_x"]
        masts["interpolated_nest"][ID]["utm_y"] = masts["raw"][ID]["utm_y"]
        masts["interpolated_nest"][ID]["i"] = ((masts["interpolated_nest"][ID]["utm_x"] -
                                                topo["interpolated_nest"]["origin_x"]) /
                                               topo["interpolated_nest"]["dx"])
        masts["interpolated_nest"][ID]["j"] = ((masts["interpolated_nest"][ID]["utm_y"] -
                                                topo["interpolated_nest"]["origin_y"]) /
                                               topo["interpolated_nest"]["dx"])

    message_done()

    # Plot the interpolated topography
    # if plotting:
    atom_plot_new(topo, "interpolated_nest", plot_number)
    plot_number += 1


# ------------------------------------------------------------------------------
#      Rotation of topography and wind turbine coordinates, if necessary
# ------------------------------------------------------------------------------
message("Rotation")

# Rotate topography
if rotation:
    message_info("Rotate topography...")
    topo["rotated"]["data"] = \
        sp.ndimage.interpolation.rotate(topo["interpolated"]["data"],
                                        angle, axes=(1, 0),
                                        reshape=True, output=None, order=3,
                                        mode='constant',
                                        cval=topo["interpolated"]["data"].min(), # noqa
                                        prefilter=False)

    # Size of rotated topography array
    topo["rotated"]["dx"] = dx
    topo["rotated"]["nx"] = topo["rotated"]["data"].shape[1]
    topo["rotated"]["ny"] = topo["rotated"]["data"].shape[0]
    topo["rotated"]["Lx"] = topo["rotated"]["nx"] * topo["rotated"]["dx"]
    topo["rotated"]["Ly"] = topo["rotated"]["ny"] * topo["rotated"]["dx"]

    message_done()

else:
    message_warning("Rotation of topography is skipped...")
    topo["rotated"] = topo["interpolated"]


message_info("Rotate turbine/mast coordinates")

# Calculate offset of rotation centers
topo["interpolated"]["center"] = np.array([topo["interpolated"]["nx"] / 2,
                                           topo["interpolated"]["ny"] / 2])
topo["rotated"]["center"] = np.array([topo["rotated"]["nx"] / 2,
                                      topo["rotated"]["ny"] / 2])

offset_rot = topo["rotated"]["center"] - topo["interpolated"]["center"]

# Rotation of turbine coordinates
for i in range(1, n_turbines + 1):
    ID = str(i)
    turbines["rotated"][ID]["i"], turbines["rotated"][ID]["j"] = \
        Rotate2D([turbines["interpolated"][ID]["i"],
                  turbines["interpolated"][ID]["j"]],
                 topo["interpolated"]["center"], math.radians(angle))

    turbines["rotated"][ID]["i"] = int(turbines["rotated"][ID]["i"] +
                                       offset_rot[0])

    turbines["rotated"][ID]["j"] = int(turbines["rotated"][ID]["j"] +
                                       offset_rot[1])

# Rotation of mast coordinates
for i in range(1, n_masts + 1):
    ID = str(i)
    masts["rotated"][ID]["i"], masts["rotated"][ID]["j"] = \
        Rotate2D([masts["interpolated"][ID]["i"],
                  masts["interpolated"][ID]["j"]],
                 topo["interpolated"]["center"], math.radians(angle))

    masts["rotated"][ID]["i"] = masts["rotated"][ID]["i"] + offset_rot[0]
    masts["rotated"][ID]["j"] = masts["rotated"][ID]["j"] + offset_rot[1]


# Get min/max coordinates of wind turbines after rotation,
# for cropping model domain

# First store the coordinates in numpy arrays...
x_values = np.zeros(n_turbines)
y_values = np.zeros(n_turbines)

for i in range(1, n_turbines + 1):
    ID = str(i)
    x_values[i - 1] = turbines["rotated"][ID]["i"]
    y_values[i - 1] = turbines["rotated"][ID]["j"]

# ... then get min/max values
wka_x_min = np.amin(x_values)
wka_x_max = np.amax(x_values)
wka_y_min = np.amin(y_values)
wka_y_max = np.amax(y_values)

message_done()


# ------------------------------------------------------------------------------
#      Rotation of topography and wind turbine coordinates (nesting)
# ------------------------------------------------------------------------------
if nesting:
    message("Rotation (nest)")

    # Rotate topography
    if rotation:
        message_info("Rotate topography (nest)...")
        topo["rotated_nest"]["data"] = \
            sp.ndimage.interpolation.rotate(topo["interpolated_nest"]["data"],
                                            angle, axes=(1, 0),
                                            reshape=True, output=None, order=3,
                                            mode='constant',
                                            cval=topo["interpolated_nest"]["data"].min(), # noqa
                                            prefilter=False)

        # Size of rotated topography array
        topo["rotated_nest"]["dx"] = dx
        topo["rotated_nest"]["nx"] = topo["rotated_nest"]["data"].shape[1]
        topo["rotated_nest"]["ny"] = topo["rotated_nest"]["data"].shape[0]
        topo["rotated_nest"]["Lx"] = topo["rotated_nest"]["nx"] * topo["rotated_nest"]["dx"]
        topo["rotated_nest"]["Ly"] = topo["rotated_nest"]["ny"] * topo["rotated_nest"]["dx"]

        message_done()

    else:
        message_warning("Rotation of topography is skipped...")
        topo["rotated_nest"] = topo["interpolated_nest"]

    message_info("Rotate turbine/mast coordinates (nest)")

    # Calculate offset of rotation centers
    topo["interpolated_nest"]["center"] = np.array([topo["interpolated_nest"]["nx"] / 2,
                                                    topo["interpolated_nest"]["ny"] / 2])
    topo["rotated_nest"]["center"] = np.array([topo["rotated_nest"]["nx"] / 2,
                                               topo["rotated_nest"]["ny"] / 2])

    offset_rot_nest = topo["rotated_nest"]["center"] - topo["interpolated_nest"]["center"]

    # Rotation of turbine coordinates
    for i in range(1, n_turbines + 1):
        ID = str(i)
        turbines["rotated_nest"][ID]["i"], turbines["rotated_nest"][ID]["j"] = \
            Rotate2D([turbines["interpolated_nest"][ID]["i"],
                      turbines["interpolated_nest"][ID]["j"]],
                     topo["interpolated_nest"]["center"], math.radians(angle))

        turbines["rotated_nest"][ID]["i"] = int(turbines["rotated_nest"][ID]["i"] +
                                                offset_rot_nest[0])

        turbines["rotated_nest"][ID]["j"] = int(turbines["rotated_nest"][ID]["j"] +
                                                offset_rot_nest[1])

    # Rotation of mast coordinates
    for i in range(1, n_masts + 1):
        ID = str(i)
        masts["rotated_nest"][ID]["i"], masts["rotated_nest"][ID]["j"] = \
            Rotate2D([masts["interpolated_nest"][ID]["i"],
                      masts["interpolated_nest"][ID]["j"]],
                     topo["interpolated_nest"]["center"], math.radians(angle))

        masts["rotated_nest"][ID]["i"] = masts["rotated_nest"][ID]["i"] + offset_rot_nest[0]
        masts["rotated_nest"][ID]["j"] = masts["rotated_nest"][ID]["j"] + offset_rot_nest[1]

    # Get min/max coordinates of wind turbines after rotation,
    # for cropping model domain

    # First store the coordinates in numpy arrays...
    x_values_nest = np.zeros(n_turbines)
    y_values_nest = np.zeros(n_turbines)

    for i in range(1, n_turbines + 1):
        ID = str(i)
        x_values_nest[i - 1] = turbines["rotated_nest"][ID]["i"]
        y_values_nest[i - 1] = turbines["rotated_nest"][ID]["j"]

    # ... then get min/max values
    wka_x_min_nest = np.amin(x_values_nest)
    wka_x_max_nest = np.amax(x_values_nest)
    wka_y_min_nest = np.amin(y_values_nest)
    wka_y_max_nest = np.amax(y_values_nest)

    message_done()


# --------------------------------
# Temporary: some hardcoded stuff... :-(
cs_distances = [50, 500]
cs_lengths = [500, 500]

# For the parent domain
for i in range(1, n_crosssections + 1):
    ID = str(i)
    # print(crosssection_1_d / topo["rotated"]["dx"])
    crosssections["rotated"][ID]["origin_i"] = (436410.0 - topo["interpolated"]["origin_x"]) / dx
    crosssections["rotated"][ID]["origin_j"] = (4743717.5 - topo["interpolated"]["origin_y"]) / dx
    crosssections["rotated"][ID]["center_i"] = turbines["rotated"]["1"]["i"] - int(cs_distances[i - 1] / topo["rotated"]["dx"])
    crosssections["rotated"][ID]["center_j"] = turbines["rotated"]["1"]["j"]
    crosssections["rotated"][ID]["north_i"] = crosssections["rotated"][ID]["center_i"]
    crosssections["rotated"][ID]["north_j"] = crosssections["rotated"][ID]["center_j"] + int(0.5 * cs_lengths[i - 1] / topo["rotated"]["dx"])
    crosssections["rotated"][ID]["south_i"] = crosssections["rotated"][ID]["center_i"]
    crosssections["rotated"][ID]["south_j"] = crosssections["rotated"][ID]["center_j"] - int(0.5 * cs_lengths[i - 1] / topo["rotated"]["dx"])

for i in range(1, n_crosssections + 1):
    ID = str(i)
    positions = ["center", "north", "south", "origin"]
    for pos in positions:
        crosssections["interpolated"][ID][pos + "_i"] = crosssections["rotated"][ID][pos + "_i"] - offset_rot[0]
        crosssections["interpolated"][ID][pos + "_j"] = crosssections["rotated"][ID][pos + "_j"] - offset_rot[1]

    positions = ["center", "north", "south", "origin"]
    for pos in positions:
        crosssections["interpolated"][ID][pos + "_i"], crosssections["interpolated"][ID][pos + "_j"] \
            = Rotate2D([crosssections["interpolated"][ID][pos + "_i"], crosssections["interpolated"][ID][pos + "_j"]], topo["interpolated"]["center"], math.radians(-angle))

    crosssections["interpolated"][ID]["north_utm_x"] = crosssections["interpolated"][ID]["north_i"] * topo["interpolated"]["dx"] + topo["interpolated"]["origin_x"]
    crosssections["interpolated"][ID]["north_utm_y"] = crosssections["interpolated"][ID]["north_j"] * topo["interpolated"]["dx"] + topo["interpolated"]["origin_y"]
    crosssections["interpolated"][ID]["south_utm_x"] = crosssections["interpolated"][ID]["south_i"] * topo["interpolated"]["dx"] + topo["interpolated"]["origin_x"]
    crosssections["interpolated"][ID]["south_utm_y"] = crosssections["interpolated"][ID]["south_j"] * topo["interpolated"]["dx"] + topo["interpolated"]["origin_y"]
    crosssections["interpolated"][ID]["origin_utm_x"] = crosssections["interpolated"][ID]["origin_i"] * topo["interpolated"]["dx"] + topo["interpolated"]["origin_x"]
    crosssections["interpolated"][ID]["origin_utm_y"] = crosssections["interpolated"][ID]["origin_j"] * topo["interpolated"]["dx"] + topo["interpolated"]["origin_y"]

# Again the same for the nest domain
if nesting:

    for i in range(1, n_crosssections + 1):
        ID = str(i)
        # print(crosssection_1_d / topo["rotated_nest"]["dx"])
        crosssections["rotated_nest"][ID]["origin_i"] = (436410.0 - topo["interpolated_nest"]["origin_x"]) / dx
        crosssections["rotated_nest"][ID]["origin_j"] = (4743717.5 - topo["interpolated_nest"]["origin_y"]) / dx
        crosssections["rotated_nest"][ID]["center_i"] = turbines["rotated_nest"]["8"]["i"] - int(cs_distances[i - 1] / topo["rotated_nest"]["dx"])
        crosssections["rotated_nest"][ID]["center_j"] = turbines["rotated_nest"]["8"]["j"]
        # print(crosssections["rotated_nest"][ID]["center_i"])
        # print(crosssections["rotated_nest"][ID]["center_j"])
        crosssections["rotated_nest"][ID]["north_i"] = crosssections["rotated_nest"][ID]["center_i"]
        crosssections["rotated_nest"][ID]["north_j"] = crosssections["rotated_nest"][ID]["center_j"] + int(0.5 * cs_lengths[i - 1] / topo["rotated_nest"]["dx"])
        crosssections["rotated_nest"][ID]["south_i"] = crosssections["rotated_nest"][ID]["center_i"]
        crosssections["rotated_nest"][ID]["south_j"] = crosssections["rotated_nest"][ID]["center_j"] - int(0.5 * cs_lengths[i - 1] / topo["rotated_nest"]["dx"])

    for i in range(1, n_crosssections + 1):
        ID = str(i)
        positions = ["center", "north", "south", "origin"]
        for pos in positions:
            crosssections["interpolated_nest"][ID][pos + "_i"] = crosssections["rotated"][ID][pos + "_i"] - offset_rot[0]
            crosssections["interpolated_nest"][ID][pos + "_j"] = crosssections["rotated"][ID][pos + "_j"] - offset_rot[1]

        positions = ["center", "north", "south", "origin"]
        for pos in positions:
            crosssections["interpolated_nest"][ID][pos + "_i"], crosssections["interpolated_nest"][ID][pos + "_j"] \
                = Rotate2D([crosssections["interpolated_nest"][ID][pos + "_i"], crosssections["interpolated_nest"][ID][pos + "_j"]], topo["interpolated_nest"]["center"], math.radians(-angle))

        crosssections["interpolated_nest"][ID]["north_utm_x"] = crosssections["interpolated_nest"][ID]["north_i"] * topo["interpolated_nest"]["dx"] + topo["interpolated_nest"]["origin_x"]
        crosssections["interpolated_nest"][ID]["north_utm_y"] = crosssections["interpolated_nest"][ID]["north_j"] * topo["interpolated_nest"]["dx"] + topo["interpolated_nest"]["origin_y"]
        crosssections["interpolated_nest"][ID]["south_utm_x"] = crosssections["interpolated_nest"][ID]["south_i"] * topo["interpolated_nest"]["dx"] + topo["interpolated_nest"]["origin_x"]
        crosssections["interpolated_nest"][ID]["south_utm_y"] = crosssections["interpolated_nest"][ID]["south_j"] * topo["interpolated_nest"]["dx"] + topo["interpolated_nest"]["origin_y"]
        crosssections["interpolated_nest"][ID]["origin_utm_x"] = crosssections["interpolated_nest"][ID]["origin_i"] * topo["interpolated_nest"]["dx"] + topo["interpolated_nest"]["origin_x"]
        crosssections["interpolated_nest"][ID]["origin_utm_y"] = crosssections["interpolated_nest"][ID]["origin_j"] * topo["interpolated_nest"]["dx"] + topo["interpolated_nest"]["origin_y"]


# # Plot the interpolated topography AGAIN! ;-)
# if plotting:
#     atom_plot_new(topo, "interpolated", plot_number)
#     plot_number += 1

# Plot the rotated topography
if plotting:
    atom_plot_new(topo, "rotated", plot_number)
    plot_number += 1

# if nesting:
#     # Plot the interpolated topography AGAIN! ;-) (nest)
#     if plotting:
#         atom_plot_new(topo, "interpolated_nest", plot_number)
#         plot_number += 1
#
#     # Plot the rotated topography (nest)
#     if plotting:
#         atom_plot_new(topo, "rotated_nest", plot_number)
#         plot_number += 1


# ------------------------------------------------------------------------------
#                  Calculate indices of model domain
# ------------------------------------------------------------------------------
message("Calculate indices of model domain")

# Convert lengths and distances to number of gridpoints
nx_inflow = lx2nx(Lx_inflow, dx)
nx_outflow = lx2nx(Lx_outflow, dx)
nx_smooth_left = lx2nx(Lx_smooth_left, dx)
nx_smooth_right = lx2nx(Lx_smooth_right, dx)
ny_cyclic = lx2nx(Ly_cyclic, dx)
nx_d_left = lx2nx(Dx_left, dx)
nx_d_right = lx2nx(Dx_right, dx)
ny_d = lx2nx(Dy, dx)

# Horizontal extensions of wind turbines (depends on wind direction)
nx_turbines = int(wka_x_max) - int(wka_x_min)
ny_turbines = int(wka_y_max) - int(wka_y_min)

Lx_turbines = nx_turbines * topo["rotated"]["dx"]

# Calculate (preliminary) model domain size for PALM,
# will be slightly modified later to fit prescribed processor topology
message("Preliminary model domain size")

nx_tmp = nx_inflow + nx_smooth_left + nx_d_left + nx_turbines \
    + nx_d_right + nx_smooth_right + nx_outflow

nx_crop = nx_smooth_left + nx_d_left + nx_turbines \
    + nx_d_right + nx_smooth_right

ny_tmp = ny_cyclic + ny_d * 2 + ny_turbines


# Do the same for the parameters of the nesting domain
if nesting:
    nx_d_left_nest = lx2nx(Dx_left_nest, dx_nest)
    nx_d_right_nest = lx2nx(Dx_right_nest, dx_nest)
    ny_d_nest = lx2nx(Dy_nest, dx_nest)

    # Horizontal extensions of wind turbines (depends on wind direction
    nx_turbines_nest = int(wka_x_max_nest) - int(wka_x_min_nest)
    ny_turbines_nest = int(wka_y_max_nest) - int(wka_y_min_nest)

    Lx_turbines_nest = nx_turbines_nest * topo["rotated_nest"]["dx"]

    nx_crop_nest = nx_d_left_nest + nx_turbines_nest + nx_d_right_nest
    ny_crop_nest = 2 * ny_d_nest + ny_turbines_nest


# ------------------------------------------------------------------------------
#                   Adjust model domain size for PALM
# ------------------------------------------------------------------------------
message_info("Processor grid:")
npex = int(nx_tmp / ngp_x + 1)
print("     npex = " + str(npex))
npey = int(ny_tmp / ngp_y + 1)
print("     npey = " + str(npey))
nx_new = npex * ngp_x
ny_new = npey * ngp_y

message_info("Adjust model domain size for nx")
nx_inflow = nx_inflow + nx_new - nx_tmp
diff_x = nx_new - nx_tmp
nx_tmp = nx_new
Lx_inflow_old = Lx_inflow
Lx_inflow = nx_inflow * dx
print("     Added " + str(diff_x) + " gridpoints in x-direction")

message_info("Adjust model domain size for ny")
ny_turbines = ny_turbines + ny_new - ny_tmp
diff_y = ny_new - ny_tmp
ny_tmp = ny_new
Ly_turbines = ny_turbines * dx
print("     Added " + str(diff_y) + " gridpoints in y-direction")

message("Final PALM model domain size")
print("     nx: " + str(nx_tmp) + ", Lx: " + str(nx_tmp * dx) + " m")
print("     ny: " + str(ny_tmp) + ", Ly: " + str(ny_tmp * dx) + " m")

procx = float(nx_tmp) / float(npex)
procy = float(ny_tmp) / float(npey)
print("     nx/npex = " + str(procx) + ", with npex = " + str(npex))
print("     ny/npey = " + str(procy) + ", with npey = " + str(npey))


# ------------------------------------------------------------------------------
#                   Adjust model domain size for PALM (nest)
# ------------------------------------------------------------------------------
if nesting:
    message("Preliminary model domain size (nest)")
    message_info("Processor grid (nest):")
    npex_nest = int(nx_crop_nest / ngp_x + 1)
    print("     npex = " + str(npex))
    npey_nest = int(ny_crop_nest / ngp_y + 1)
    print("     npey = " + str(npey))
    nx_new_nest = npex_nest * ngp_x
    ny_new_nest = npey_nest * ngp_y

    message_info("Adjust model domain size for nx (nest)")
    print("  old nx_crop_nest: " + str(nx_crop_nest))
    diff_x_nest = nx_new_nest - nx_crop_nest
    nx_crop_nest = nx_new_nest
    print("  new nx_crop_nest: " + str(nx_crop_nest))
    print("     Added " + str(diff_x_nest) + " gridpoints in x-direction (nest)")

    message_info("Adjust model domain size for ny (nest)")
    print("  old ny_crop_nest: " + str(ny_crop_nest))
    diff_y_nest = ny_new_nest - ny_crop_nest
    ny_crop_nest = ny_new_nest
    print("  new ny_crop_nest: " + str(ny_crop_nest))
    print("     Added " + str(diff_y_nest) + " gridpoints in y-direction (nest)")

    message("Final PALM model domain size (nest)")
    print("     nx: " + str(nx_crop_nest) + ", Lx: " + str(nx_crop_nest * dx_nest) + " m")
    print("     ny: " + str(ny_crop_nest) + ", Ly: " + str(ny_crop_nest * dx_nest) + " m")

    procx_nest = float(nx_crop_nest) / float(npex_nest)
    procy_nest = float(ny_crop_nest) / float(npey_nest)
    print("     nx/npex = " + str(procx_nest) + ", with npex = " + str(npex_nest))
    print("     ny/npey = " + str(procy_nest) + ", with npey = " + str(npey_nest))


# ------------------------------------------------------------------------------
#                Coordinates for cutting out model domain
# ------------------------------------------------------------------------------
x0_cut = int(wka_x_min) - nx_d_left - nx_smooth_left  # - nx_inflow
y0_cut = int(wka_y_min) - ny_d - ny_cyclic
xmax_cut = x0_cut + nx_crop
ymax_cut = y0_cut + ny_tmp  # + ny_cyclic

for i in range(1, n_turbines + 1):
    ID = str(i)
    turbines["cropped"][ID]["i"] = turbines["rotated"][ID]["i"] - x0_cut
    turbines["cropped"][ID]["j"] = turbines["rotated"][ID]["j"] - y0_cut
    turbines["PALM"][ID]["i"] = turbines["rotated"][ID]["i"] - x0_cut + nx_inflow
    turbines["PALM"][ID]["j"] = turbines["rotated"][ID]["j"] - y0_cut

for i in range(1, n_masts + 1):
    ID = str(i)
    masts["cropped"][ID]["i"] = masts["rotated"][ID]["i"] - x0_cut
    masts["cropped"][ID]["j"] = masts["rotated"][ID]["j"] - y0_cut
    masts["PALM"][ID]["i"] = masts["rotated"][ID]["i"] - x0_cut + nx_inflow
    masts["PALM"][ID]["j"] = masts["rotated"][ID]["j"] - y0_cut

for i in range(1, n_crosssections + 1):
    ID = str(i)
    positions = ["center", "north", "south", "origin"]
    for pos in positions:
        crosssections["cropped"][ID][pos + "_i"] = crosssections["rotated"][ID][pos + "_i"] - x0_cut
        crosssections["cropped"][ID][pos + "_j"] = crosssections["rotated"][ID][pos + "_j"] - y0_cut
        crosssections["PALM"][ID][pos + "_i"] = crosssections["rotated"][ID][pos + "_i"] - x0_cut + nx_inflow
        crosssections["PALM"][ID][pos + "_j"] = crosssections["rotated"][ID][pos + "_j"] - y0_cut


# ------------------------------------------------------------------------------
#                Coordinates for cutting out model domain (nest)
# ------------------------------------------------------------------------------
if nesting:
    x0_cut_nest = int(wka_x_min_nest) - nx_d_left_nest
    y0_cut_nest = int(wka_y_min_nest) - ny_d_nest
    xmax_cut_nest = x0_cut_nest + nx_crop_nest
    ymax_cut_nest = y0_cut_nest + ny_crop_nest

    for i in range(1, n_turbines + 1):
        ID = str(i)
        turbines["cropped_nest"][ID]["i"] = turbines["rotated_nest"][ID]["i"] - x0_cut_nest
        turbines["cropped_nest"][ID]["j"] = turbines["rotated_nest"][ID]["j"] - y0_cut_nest
        turbines["PALM_nest"][ID]["i"] = turbines["rotated_nest"][ID]["i"] - x0_cut_nest
        turbines["PALM_nest"][ID]["j"] = turbines["rotated_nest"][ID]["j"] - y0_cut_nest

    for i in range(1, n_masts + 1):
        ID = str(i)
        masts["cropped_nest"][ID]["i"] = masts["rotated_nest"][ID]["i"] - x0_cut_nest
        masts["cropped_nest"][ID]["j"] = masts["rotated_nest"][ID]["j"] - y0_cut_nest
        masts["PALM_nest"][ID]["i"] = masts["rotated_nest"][ID]["i"] - x0_cut_nest
        masts["PALM_nest"][ID]["j"] = masts["rotated_nest"][ID]["j"] - y0_cut_nest


# ------------------------------------------------------------------------------
#                      Difference between parent and nest
# ------------------------------------------------------------------------------
if nesting:
    x_diff = (nx_inflow + nx_smooth_left + nx_d_left) * dx - nx_d_left_nest * dx_nest
    y_diff = (ny_cyclic + ny_d) * dx - ny_d_nest * dx_nest

    print(" xdiff fuer nesting: " + str(x_diff))
    print(" ydiff fuer nesting: " + str(y_diff))

    crosssections["PALM"]["1"]["nest_origin_i"] = x_diff / dx
    crosssections["PALM"]["1"]["nest_origin_j"] = y_diff / dx

# ------------------------------------------------------------------------------
#                          Cut out topography data
# ------------------------------------------------------------------------------
message("Cutting out topography data")

message_info("Coordinate of lower left corner of input topography data")
print("     x: " + str(zz[0]))
print("     y: " + str(zz[1]))

# Cut out topography
message_info("Cut out")

topo["cropped"]["data"] = topo["rotated"]["data"][y0_cut:ymax_cut, x0_cut:xmax_cut]
topo["cropped"]["dx"] = dx
topo["cropped"]["nx"] = topo["cropped"]["data"].shape[1]
topo["cropped"]["ny"] = topo["cropped"]["data"].shape[0]
topo["cropped"]["Lx"] = topo["cropped"]["nx"] * topo["cropped"]["dx"]
topo["cropped"]["Ly"] = topo["cropped"]["ny"] * topo["cropped"]["dx"]

message_done()

# Plot the cropped topography
if plotting:
    atom_plot_new(topo, "cropped", plot_number)
    plot_number += 1

topo["PALM"]["origin_x"] = zz[0] + (x0_cut * dx) - (nx_inflow * dx)
topo["PALM"]["origin_y"] = zz[1] + (y0_cut * dx)


# ------------------------------------------------------------------------------
#                          Cut out topography data
# ------------------------------------------------------------------------------
if nesting:
    message("Cutting out topography data (nest)")

    message_info("Coordinate of lower left corner of input topography data")
    print("     x: " + str(zz[0]))
    print("     y: " + str(zz[1]))

    # Cut out topography
    message_info("Cut out (nest)")

    topo["cropped_nest"]["data"] = topo["rotated_nest"]["data"][y0_cut_nest:ymax_cut_nest, x0_cut_nest:xmax_cut_nest]
    topo["cropped_nest"]["dx"] = dx_nest
    topo["cropped_nest"]["nx"] = topo["cropped_nest"]["data"].shape[1]
    topo["cropped_nest"]["ny"] = topo["cropped_nest"]["data"].shape[0]
    topo["cropped_nest"]["Lx"] = topo["cropped_nest"]["nx"] * topo["cropped_nest"]["dx"]
    topo["cropped_nest"]["Ly"] = topo["cropped_nest"]["ny"] * topo["cropped_nest"]["dx"]

    message_done()

    # Plot the cropped topography
    if plotting:
        atom_plot_new(topo, "cropped_nest", plot_number)
        plot_number += 1

    # possible TODO: replace zz here???
    topo["PALM_nest"]["origin_x"] = zz[0] + (x0_cut_nest * dx_nest)
    topo["PALM_nest"]["origin_y"] = zz[1] + (y0_cut_nest * dx_nest)


# ------------------------------------------------------------------------------
#                 Coordinates of crosssection yz_crosssections
# ------------------------------------------------------------------------------
if yz_crosssections:
    zz_cs_1 = np.array([0, 0])
    zz_cs_1[0] = topo["PALM"]["origin_x"] + (crosssection_1 * dx)
    zz_cs_1[1] = topo["PALM"]["origin_y"] + ((crosssection_y - crosssection_1_d / 2) * dx)

    zz_cs_2 = np.array([0, 0])
    zz_cs_2[0] = topo["PALM"]["origin_x"] + (crosssection_2 * dx)
    zz_cs_2[1] = topo["PALM"]["origin_y"] + ((crosssection_y - crosssection_2_d / 2) * dx)

    message_info("ZZ PALM 1,2")
    print(topo["PALM"]["origin_x"])
    print(topo["PALM"]["origin_y"])

    message_info("ZZ CS 1")
    print(zz_cs_1[0])
    print(zz_cs_1[1])

    message_info("ZZ CS 2")
    print(zz_cs_2[0])
    print(zz_cs_2[1])

    message_info("SECURITY CHECK CS1")
    print(zz_cs_1[0] - topo["PALM"]["origin_x"])
    print(zz_cs_1[1] - topo["PALM"]["origin_y"])

    message_info("SECURITY CHECK CS2")
    print(zz_cs_2[0] - topo["PALM"]["origin_x"])
    print(zz_cs_2[1] - topo["PALM"]["origin_y"])


# ------------------------------------------------------------------------------
#                        Modify topography data
# ------------------------------------------------------------------------------
message("Modify topography data")

# This seems to be UTM, however we are in a rotated system now!!! :-(
topo["PALM"]["origin_x"] = zz[0] + (x0_cut * dx) - (nx_inflow * dx)
topo["PALM"]["origin_y"] = zz[1] + (y0_cut * dx)

topo["PALM"]["nx"] = nx_tmp
topo["PALM"]["ny"] = topo["cropped"]["ny"]
topo["PALM"]["dx"] = dx
topo["PALM"]["Lx"] = topo["PALM"]["nx"] * topo["PALM"]["dx"]
topo["PALM"]["Ly"] = topo["PALM"]["ny"] * topo["PALM"]["dx"]

# Initialize final data array
topo["PALM"]["data"] = np.zeros((topo["PALM"]["ny"], topo["PALM"]["nx"]))
# Fill array with cropped topography
topo["PALM"]["data"][:, nx_inflow:nx_inflow + topo["cropped"]["nx"]] = topo["cropped"]["data"][:, :]

# Left transition subdomain (linear blending from flat to unmodified topography)
message_info("Crossfading left transition")
for i in range(nx_inflow, nx_inflow + nx_smooth_left):
    sum_tmp = 0
    fac = 1.0 * (i - nx_inflow) / nx_smooth_left
    for j in range(0, ny_tmp):
        sum_tmp = sum_tmp + topo["PALM"]["data"][j, i]
    sum_tmp = sum_tmp / ny_tmp
    topo["PALM"]["data"][:, i] = fac * topo["PALM"]["data"][:, i] + (1.0 - fac) * sum_tmp
    if i == nx_inflow:
        inflow_height = sum_tmp
message_done()

# Inflow plane
message_info("Inflow plane")
for i in range(0, nx_inflow):
    topo["PALM"]["data"][:, i] = inflow_height
message_done()

# Right transition subdomain (linear blending from unmodified topography to flat)
message_info("Crossfading right transition")
for i in range(nx_inflow + nx_crop - nx_smooth_right, nx_inflow + nx_crop):
    sum_tmp = 0
    fac = 1.0 * (i - (nx_inflow + nx_crop - nx_smooth_right)) / (nx_smooth_right)
    for j in range(0, ny_tmp):
        sum_tmp = sum_tmp + topo["PALM"]["data"][j, i]
    sum_tmp = sum_tmp / ny_tmp
    topo["PALM"]["data"][:, i] = (1.0 - fac) * topo["PALM"]["data"][:, i] + fac * sum_tmp
    if i == nx_inflow + nx_crop - 1:
        outflow_height = sum_tmp
message_done()

# Outflow plane
message_info("Outflow plane")
for i in range(nx_inflow + nx_crop, nx_inflow + nx_crop + nx_outflow):
    topo["PALM"]["data"][:, i] = outflow_height
message_done()

# Plot the almost final topography
if plotting:
    atom_plot_new(topo, "PALM", plot_number)
    plot_number += 1

# Make periodic in y-direction
message_info("Make topography periodic in y-direction")
for i in range(nx_inflow, nx_inflow + nx_crop):
    for j in range(0, ny_cyclic):
        fac = float(j) / ny_cyclic
        topo["PALM"]["data"][j, i] = (1.0 - fac) * topo["PALM"]["data"][j + ny_tmp - ny_cyclic, i] + fac * topo["PALM"]["data"][j, i]
message_done()

# Plot the topography for usage in PALM, lowest point NOT set to zero yet
if plotting:
    atom_plot_new(topo, "PALM", plot_number)
    plot_number += 1

# Set lowest point to zero
message_info("Set lowest point to zero")
min_value = topo["PALM"]["data"].min()
topo["PALM"]["data"] = topo["PALM"]["data"] - min_value
message_done()

if nesting:
    topo["PALM_nest"] = topo["cropped_nest"]
    topo["PALM_nest"]["data"] = topo["PALM_nest"]["data"] - min_value
    # print(topo["PALM_nest"]["data"])
    # print(topo["PALM"]["data"])

# Plot the final topography for usage in PALM
atom_plot_new(topo, "PALM", plot_number)
if nesting:
    atom_plot_new(topo, "PALM_nest", plot_number)

# Update heights of inflow and outflow planes
inflow_height = inflow_height - min_value
outflow_height = outflow_height - min_value
print("  Minimum value: " + str(min_value))
print("  Inflow height: " + str(inflow_height))
print("  Outflow height: " + str(outflow_height))

# generate 2 2d grids for the x & y bounds
topo["PALM"]["dx"] = dx
topo["PALM"]["Lx"] = topo["PALM"]["nx"] * topo["PALM"]["dx"]
topo["PALM"]["Ly"] = topo["PALM"]["ny"] * topo["PALM"]["dx"]

message("Wind turbine coordinates for usage in PALM:")
print("replace old message")
# wka_ij_palm = np.around(wka_ij_rot, decimals=0)
# wka_ij_palm_int = wka_ij_palm.astype(int)

# for turbine in range(n_turbines):
print("replace old message")
# print("  Turbine #" + str(turbine+1) + " | x: " + str(wka_ij_palm[turbine,0]) + ", y: " + str(wka_ij_palm[turbine,1]))

message("Mast coordinates for usage in PALM:")
# masts_ij_palm = np.around(masts_ij_rot, decimals=0)
# masts_ij_palm_int = masts_ij_palm.astype(int)

# for mast in range(n_masts):
print("replace old message")
# print("  Mast #" + str(mast+1) + " | x: " + str(masts_ij_palm[mast,0]) + ", y: " + str(masts_ij_palm[mast,1]))


# ------------------------------------------------------------------------------
#             Calculate indices of yz_crosssections
#
#                       -->  clean this section later
# ------------------------------------------------------------------------------
if yz_crosssections:

    message("yz_crosssections")
    print(zz)

    gp_x0_1 = crosssection_1
    gp_y0_1 = crosssection_y - crosssection_1_d / 2
    gp_x_1 = 0
    gp_y_1 = crosssection_1_d

    message_info("Indices of CS 1")
    print('x ' + str(gp_x0_1))
    print('y0 ' + str(gp_y0_1))
    print('ye ' + str(gp_y0_1 + gp_y_1))

    gp_x0_2 = crosssection_2
    gp_y0_2 = crosssection_y - crosssection_2_d / 2
    gp_x_2 = 0
    gp_y_2 = crosssection_2_d

    message_info("Indices of CS 2")
    print('x ' + str(gp_x0_2))
    print('y0 ' + str(gp_y0_2))
    print('ye ' + str(gp_y0_2 + gp_y_2))


# ------------------------------------------------------------------------------
#                         Output new ASCII file
# ------------------------------------------------------------------------------
message("Output new ASCII file(s)...")

outputfile = output_dir + 'INPUT/' + case_name + "_topo"
print(outputfile)

# Flip data upside down
topo["PALM"]["data"] = np.flipud(topo["PALM"]["data"])
np.savetxt(outputfile, topo["PALM"]["data"], fmt='%1.3f')


# ------------------------------------------------------------------------------
#                         Output new ASCII file (nest)
# ------------------------------------------------------------------------------
if nesting:
    outputfile = output_dir + 'INPUT/' + case_name + "_02_topo"
    print(outputfile)

    # Flip data upside down
    topo["PALM_nest"]["data"] = np.flipud(topo["PALM_nest"]["data"])
    np.savetxt(outputfile, topo["PALM_nest"]["data"], fmt='%1.3f')


# ------------------------------------------------------------------------------
#                         Create the p3d-File
# ------------------------------------------------------------------------------

if not nesting:
    # Path to p3d-file which will be modified at the end of this program
    p3d_path = output_dir + 'INPUT/' + case_name + "_p3d"

    # Copy p3d-template to output-directory
    copyfile(p3d_template, p3d_path)

    nxm = nx_tmp - 1
    nym = ny_tmp - 1
    nz_dict = {40: 84, 20: 168, 10: 336, 5: 672, 2.5: 9999}
    nzm = nz_dict[dx]

    dy = dx
    dz = dx

    x_values = np.zeros(n_turbines, dtype=np.int)
    y_values = np.zeros(n_turbines, dtype=np.int)

    for i in range(1, n_turbines + 1):
        ID = str(i)
        x_values[i - 1] = int(turbines["PALM"][ID]["i"])
        y_values[i - 1] = int(turbines["PALM"][ID]["j"])

    wka_x_str = ",".join(map(str, x_values))
    wka_y_str = ",".join(map(str, y_values))

    p3d_var = np.array([nxm, nym, nzm, dx, dy, dz, npex, npey, n_turbines, wka_x_str, wka_y_str])
    p3d_varnames = np.array(["nx", "ny", "nz", "dx", "dy", "dz", "npex", "npey", "n_turbines", "wka_x", "wka_y"])

    for i in range(len(p3d_var)):
        replace(p3d_path, "<" + p3d_varnames[i] + ">", str(p3d_var[i]))

else:
    # Path to p3d-file which will be modified at the end of this program
    p3d_parent_path = output_dir + 'INPUT/' + case_name + "_p3d"
    p3d_nest_path = output_dir + 'INPUT/' + case_name + "_02_p3d"

    # Copy p3d-template to output-directory
    copyfile(p3d_template_parent, p3d_parent_path)
    copyfile(p3d_template_nest, p3d_nest_path)

    nxm = nx_tmp - 1
    nym = ny_tmp - 1
    nz_dict = {40: 84, 20: 168, 10: 336, 5: 672, 2.5: 9999}
    nzm = nz_dict[dx]
    dy = dx
    dz = dx

    nxm_nest = nx_crop_nest - 1
    nym_nest = ny_crop_nest - 1
    nz_nest_dict = {40: 42, 20: 84, 10: 168, 5: 336, 2.5: 672}
    nzm_nest = nz_nest_dict[dx]
    dy_nest = dx_nest
    dz_nest = dx_nest

    x_values = np.zeros(n_turbines, dtype=np.int)
    y_values = np.zeros(n_turbines, dtype=np.int)

    x_values_nest = np.zeros(n_turbines, dtype=np.int)
    y_values_nest = np.zeros(n_turbines, dtype=np.int)

    for i in range(1, n_turbines + 1):
        ID = str(i)
        x_values[i - 1] = int(turbines["PALM"][ID]["i"])
        y_values[i - 1] = int(turbines["PALM"][ID]["j"])
        x_values_nest[i - 1] = int(turbines["PALM_nest"][ID]["i"])
        y_values_nest[i - 1] = int(turbines["PALM_nest"][ID]["j"])

    wka_x_str = ",".join(map(str, x_values))
    wka_y_str = ",".join(map(str, y_values))
    wka_x_nest_str = ",".join(map(str, x_values_nest))
    wka_y_nest_str = ",".join(map(str, y_values_nest))

    np_parent = npex * npey
    np_nest = npex_nest * npey_nest

    # Vars to be written in p3d file
    p3d_var = np.array([nxm, nym, nzm, dx, dy, dz, npex, npey, n_turbines, wka_x_str, wka_y_str, np_parent, np_nest, x_diff, y_diff])
    p3d_nest_var = np.array([nxm_nest, nym_nest, nzm_nest, dx_nest, dy_nest, dz_nest, npex_nest, npey_nest, n_turbines, wka_x_nest_str, wka_y_nest_str])

    # Placeholders in p3d template
    p3d_varnames = np.array(["nx", "ny", "nz", "dx", "dy", "dz", "npex", "npey",
                             "n_turbines", "wka_x", "wka_y",
                             "np_parent", "np_nest", "xdiff", "ydiff"])
    p3d_nest_varnames = np.array(["nx", "ny", "nz", "dx", "dy", "dz", "npex", "npey",
                                  "n_turbines", "wka_x", "wka_y"])

    for i in range(len(p3d_var)):
        replace(p3d_parent_path, "<" + p3d_varnames[i] + ">", str(p3d_var[i]))

    for i in range(len(p3d_nest_var)):
        replace(p3d_nest_path, "<" + p3d_varnames[i] + ">", str(p3d_nest_var[i]))

# ------------------------------------------------------------------------------
#                           Output INFO file
# ------------------------------------------------------------------------------
info_file = open(output_dir + case_name + '_info.txt', 'w')
info_file.write('------------------------------------------------------' + '\n')
info_file.write('Run-ID: ' + case_name + '\n')
info_file.write('------------------------------------------------------' + '\n')
info_file.write('dx: ' + str(dx) + ' m' + '\n')
info_file.write('Dx_left: ' + str(Dx_left) + ' m' + '\n')
info_file.write('Dx_right: ' + str(Dx_right) + ' m' + '\n')
info_file.write('Dy: ' + str(Dy) + ' m' + '\n')
info_file.write('Lx_inflow: ' + str(Lx_inflow) + ' m' + '\n')
info_file.write('Lx_inflow_old: ' + str(Lx_inflow_old) + ' m' + '\n')
info_file.write('Lx_outflow: ' + str(Lx_outflow) + ' m' + '\n')
info_file.write('Lx_smooth_left: ' + str(Lx_smooth_left) + ' m' + '\n')
info_file.write('Lx_smooth_right: ' + str(Lx_smooth_right) + ' m' + '\n')
info_file.write('Lx_turbines: ' + str(Lx_turbines) + ' m' + '\n')
info_file.write('Ly_cyclic: ' + str(Ly_cyclic) + ' m' + '\n')
info_file.write('Ly_turbines: ' + str(Ly_turbines) + ' m' + '\n')
info_file.write('------------------------------------------------------' + '\n')
info_file.write('Minimum height: ' + str(min_value) + ' m' + '\n')
info_file.write('Inflow height: ' + str(inflow_height) + ' m' + '\n')
info_file.write('Outflow height: ' + str(outflow_height) + ' m' + '\n')
info_file.write('------------------------------------------------------' + '\n')
info_file.write("Final PALM model domain size" + '\n')
info_file.write("nx: " + str(nx_tmp) + ", Lx: " + str(nx_tmp * dx) + " m" + '\n')
info_file.write("ny: " + str(ny_tmp) + ", Ly: " + str(ny_tmp * dx) + " m" + '\n')
info_file.write('------------------------------------------------------' + '\n')
info_file.write("Numerical grid:" + '\n')
info_file.write('ngp_x: ' + str(ngp_x) + '\n')
info_file.write('ngp_y: ' + str(ngp_y) + '\n')
info_file.write("nx/npex = " + str(procx) + ", with npex = " + str(npex) + '\n')
info_file.write("ny/npey = " + str(procy) + ", with npey = " + str(npey) + '\n')
info_file.write('------------------------------------------------------' + '\n')
info_file.write('--For &userpar ---------------------------------------' + '\n')
info_file.write('wka_x = ')
for turbine in range(n_turbines):
    # info_file.write(str(int(wka_ij_palm[turbine,0])) + ', ')
    info_file.write("replace with new message")
info_file.write('\n')
info_file.write('wka_y = ')
for turbine in range(n_turbines):
    # info_file.write(str(int(wka_ij_palm[turbine,1])) + ', ')
    info_file.write("replace with new message")
info_file.write('\n')
info_file.write('mm_x = ')
for mast in range(n_masts):
    # info_file.write(str(int(masts_ij_palm[mast,0])) + ', ')
    info_file.write("replace with new message")
info_file.write('\n')
info_file.write('mm_y = ')
for mast in range(n_masts):
    # info_file.write(str(int(masts_ij_palm[mast,1])) + ', ')
    info_file.write("replace with new message")
info_file.write('\n')
info_file.write('------------------------------------------------------' + '\n')

info_file.close()


# ------------------------------------------------------------------------------
#                            End of program
# ------------------------------------------------------------------------------
message_ok("ATOM finished.")

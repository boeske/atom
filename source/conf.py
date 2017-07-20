"""
--------------------------------------------------------------------------------
ATOM: Settings
Version 2.0
--------------------------------------------------------------------------------
"""

import numpy as np
import configparser
from source.messages import *
import os.path
import sys

# ------------------------------------------------------------------------------
#                               ConfigParser
# ------------------------------------------------------------------------------

config = configparser.ConfigParser()
if os.path.exists("atom.ini"):
    config.read("atom.ini")
else:
    message_error("Error: atom.ini not found!")

# Get case name and set output directory
case_name = config.get("Name", "case_name")
output_dir = "output/" + case_name + "/"

# First get number of wind turbines, then create an array to store
# coordinates of turbines
n_turbines = config.getint("WindTurbines", "n_turbines")
wka_ij = np.zeros((n_turbines, 2))
for turbine in range(n_turbines):
    coordinate_str = "coordinates_" + str(turbine + 1)
    wka_ij[turbine, :] = eval(config.get("WindTurbines", coordinate_str), {}, {})


n_masts = config.getint("Masts", "n_masts")
masts_ij = np.zeros((n_masts, 2))
for mast in range(n_masts):
    coordinate_str = "coordinates_" + str(mast + 1)
    masts_ij[mast, :] = eval(config.get("Masts", coordinate_str), {}, {})

# Plotting config
plotting = config.getboolean("Plotting", "plotting")
dpi_value = config.getint("Plotting", "dpi_value")
plot_number = 1

# Input data
inputfile = config.get("InputData", "inputfile")
dx_raw = config.getfloat("InputData", "dx_raw")
zz = np.array([config.getfloat("InputData", "origin_x"),
               config.getfloat("InputData", "origin_y")])
origin_x = config.getfloat("InputData", "origin_x")
origin_y = config.getfloat("InputData", "origin_y")
p3d_template = config.get("InputData", "p3d_template")

# Meteorology
wdir = config.getfloat("Meteorology", "wdir")

# Calculate topography rotation angle from wind direction
if wdir == 270.:
    rotation = False
    angle = 0
else:
    rotation = True
    angle = 270. - wdir

# Model domain
dx = config.getfloat("ModelDomain", "dx")
Dx_left = config.getfloat("ModelDomain", "Dx_left")
Dx_right = config.getfloat("ModelDomain", "Dx_right")
Dy = config.getfloat("ModelDomain", "Dy")
Lx_inflow = config.getfloat("ModelDomain", "Lx_inflow")
Lx_outflow = config.getfloat("ModelDomain", "Lx_outflow")
Lx_smooth_left = config.getfloat("ModelDomain", "Lx_smooth_left")
Lx_smooth_right = config.getfloat("ModelDomain", "Lx_smooth_right")
Ly_cyclic = config.getfloat("ModelDomain", "Ly_cyclic")

# Processor grid
ngp_x = config.getint("ProcessorGrid", "ngp_x")
ngp_y = config.getint("ProcessorGrid", "ngp_y")

# Nesting
nesting = config.getboolean("Nesting", "nesting")
if nesting:
    dx_nest = config.getfloat("Nesting", "dx_nest")
    p3d_template_nest = config.get("Nesting", "p3d_template_nest")
    p3d_template_parent = config.get("Nesting", "p3d_template_parent")

# Temporary defaults for nest domain size
Dx_left_nest = Dx_left / 2
Dy_nest = Dy / 2
Dx_right_nest = Dx_right / 2


# ------------------------------------------------------------------------------
#                     More variables
# ------------------------------------------------------------------------------

# Extension of wind turbines in x-dir
# (Depends on spatial distribution of wind turbines and wind direction 'wdir')
# Value will be calculated by ATOM
Lx_turbines = None

# Extension of wind turbines in y-dir
# (Depends on 'npey' and spatial distribution of wind turbines and wind
# direction 'wdir')
# Value will be calculated by ATOM
Ly_turbines = None


# ------------------------------------------------------------------------------
#                  Special settings for yz_crosssections
# ------------------------------------------------------------------------------

# Enable Output for yz_crosssections
yz_crosssections = False

crosssection_1 = 3795
crosssection_1_d = 200    # Dy = 500m
crosssection_2 = 3975
crosssection_2_d = 200    # Dy = 500m
crosssection_y = 917

n_crosssections = 2


# ------------------------------------------------------------------------------
#                  Initialize dictionarys to store all data
# ------------------------------------------------------------------------------


# Define stages of modification
stages = ["raw", "interpolated", "rotated", "cropped", "PALM"]

# Add stages for nesting if needed
if nesting:
    stages.append("interpolated_nest")
    stages.append("rotated_nest")
    stages.append("cropped_nest")
    stages.append("PALM_nest")

# This dictionary will store the topography data for each stage,
# and also metadata like the coordinate of lower left corner, and so on...
topo = {}
for stage in stages:
    topo[stage] = {}

# Keys for the turbine & mast dictionaries
keys = ["i", "j", "utm_x", "utm_y"]

# This dictionary will store information for each turbine at each mod stage
# For example coordinates, ij_coordinate, height, and so on...
turbines = {}
for stage in stages:
    turbines[stage] = {}
    for i in range(1, n_turbines + 1):
        ID = str(i)
        turbines[stage][ID] = {}
        for key in keys:
            turbines[stage][ID][key] = 0

# This dictionary will store information for each turbine at each mod stage
# For example coordinates, ij_coordinate, height, and so on...
masts = {}
for stage in stages:
    masts[stage] = {}
    for i in range(1, n_masts + 1):
        ID = str(i)
        masts[stage][ID] = {}
        for key in keys:
            masts[stage][ID][key] = 0

# This dictionary will store information for each crosssection for all(?) stages
# However, this has to be calculated backwards, beginning at stage "PALM" to "raw"!
crosssections = {}
cs_keys = ["center_i", "center_j", "north_i", "north_j", "south_i", "south_j",
           "north_utm_x", "north_utm_y", "south_utm_x", "south_utm_y",
           "origin_i", "origin_j", "origin_utm_x", "origin_utm_y",
           "nest_origin_i", "nest_origin_j", "nest_origin_utm_x", "nest_origin_utm_y"]
for stage in stages:
    crosssections[stage] = {}
    for i in range(1, n_crosssections + 1):
        ID = str(i)
        crosssections[stage][ID] = {}
        for key in cs_keys:
            crosssections[stage][ID][key] = 0

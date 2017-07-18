"""
ATOM: Functions
Version 2.0
"""

import scipy as sp
from tempfile import mkstemp
from shutil import move
from os import remove, close

# ------------------------------------------------------------------------------
#                           Definition of functions
# ------------------------------------------------------------------------------


# Function to rotate coordinates
def Rotate2D(pts, cnt, ang):
    return sp.dot(pts - cnt, sp.array([[sp.cos(-ang), sp.sin(-ang)],
                                       [-sp.sin(-ang), sp.cos(-ang)]])) + cnt


# Conversion of a length to number of gridpoints
def lx2nx(length, dx):
    nx = length / dx
    if nx < 1:
        nx = 1
    return int(nx)


# Replace 'pattern' by 'subst' in a given file
def replace(file_path, pattern, subst):
    # Create temp file
    fh, abs_path = mkstemp()
    with open(abs_path, 'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))
    close(fh)

    # Remove original file
    remove(file_path)

    # Move new file to final destination
    move(abs_path, file_path)

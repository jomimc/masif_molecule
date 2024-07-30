import os
import numpy
from subprocess import Popen, PIPE
import pymesh

from default_config.global_vars import apbs_bin, pdb2pqr_bin, multivalue_bin
import random

"""
computeAPBS.py: Wrapper function to compute the Poisson Boltzmann electrostatics for a surface using APBS.
Pablo Gainza - LPDI STI EPFL 2019
This file is part of MaSIF.
Released under an Apache License 2.0
"""

def computeAPBS(vertices, path_input, tmp_dir):
    """
        Calls APBS, pdb2pqr, and multivalue and returns the charges per vertex
    """
#   fields = tmp_file_base.split("/")[0:-1]
#   directory = "/".join(fields) + "/"
    filename_base = path_input.stem
    pdbname = path_input.name
    args = [
        pdb2pqr_bin,
        "--ff=parse",
        "--whitespace",
        "--noopt",
        "--apbs-input",
        pdbname,
        filename_base,
    ]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=str(tmp_dir))
    stdout, stderr = p2.communicate()

    args = [apbs_bin, filename_base + ".in"]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=str(tmp_dir))
    stdout, stderr = p2.communicate()

    with open(tmp_dir.joinpath(f"{filename_base}.csv"), "w") as vertfile:
        for vert in vertices:
            vertfile.write("{},{},{}\n".format(vert[0], vert[1], vert[2]))

    args = [
        multivalue_bin,
        filename_base + ".csv",
        filename_base + ".dx",
        filename_base + "_out.csv",
    ]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=str(tmp_dir))
    stdout, stderr = p2.communicate()

    # Read the charge file
    with open(tmp_dir.joinpath(f"{filename_base}_out.csv"), "r") as chargefile:
        charges = numpy.array([0.0] * len(vertices))
        for ix, line in enumerate(chargefile.readlines()):
            charges[ix] = float(line.split(",")[3])

    return charges


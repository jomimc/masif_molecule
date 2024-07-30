import numpy as np
from subprocess import Popen, PIPE

from default_config.global_vars import apbs_bin, pdb2pqr_bin, multivalue_bin

"""
computeAPBS.py: Wrapper function to compute the Poisson Boltzmann electrostatics for a surface using APBS.
Pablo Gainza - LPDI STI EPFL 2019
This file is part of MaSIF.
Released under an Apache License 2.0
"""


def convert_apbs_input(path_input, tmp_dir):
    path_template = 'apbs_input.in'
    path_apbs_input = tmp_dir.joinpath(f"{path_input.stem}.in")
    path_dx = tmp_dir.joinpath(f"{path_input.stem}")
    lines = open(path_template, 'r').readlines()
    lines[1] = f"    mol pqr {str(path_input)}\n"
    lines[23] = f"    write pot dx {str(path_dx)}\n"
    with open(path_apbs_input, 'w') as o:
        for l in lines:
            o.write(l)
        

def computeAPBS(vertices, path_input, tmp_dir):
    """
        Calls APBS, pdb2pqr, and multivalue and returns the charges per vertex
    """

    filename_base = path_input.stem
    pdbname = path_input.name

    if path_input.suffix == '.pdb':
        # Convert to PQR
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

    elif path_input.suffix == '.pqr':
        # Create APBS input
        convert_apbs_input(path_input, tmp_dir)

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
        charges = np.array([0.0] * len(vertices))
        for ix, line in enumerate(chargefile.readlines()):
            charges[ix] = float(line.split(",")[3])

    return charges


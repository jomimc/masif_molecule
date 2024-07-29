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


def rewrite_apbs_input(path_pqr, model):
    path_in = "apbs_input.in"
    path_out = PATH_BASE.joinpath("tempfiles", f"{path_pqr.stem}_{model:03d}.in")
    lines = [l for l in open(path_in, 'r')]
    lines[1] = f"    mol pqr {str(path_pqr)}\n"
    lines[23] = f"    write pot dx {str(path_pqr)}\n"
    with open(path_out, 'w') as o:
        for l in lines:
            o.write(l)
    return path_out



def computeAPBS(vertices, path_input, tmp_file_base):
    """
        Calls APBS, pdb2pqr, and multivalue and returns the charges per vertex
    """

    directory = str(PATH_BASE.joinpath("tempfiles"))
    path_apbs = rewrite_apbs_input(path_pqr, model)
    args = [apbs_bin, path_apbs]


    fields = tmp_file_base.split("/")[0:-1]
    directory = "/".join(fields) + "/"
    filename_base = tmp_file_base.split("/")[-1]
    pdbname = pdb_file.split("/")[-1]
    args = [
        pdb2pqr_bin,
        "--ff=parse",
        "--whitespace",
        "--noopt",
        "--apbs-input",
        pdbname,
        filename_base,
    ]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=directory)
    stdout, stderr = p2.communicate()

    args = [apbs_bin, filename_base + ".in"]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=directory)
    stdout, stderr = p2.communicate()

    vertfile = open(directory + "/" + filename_base + ".csv", "w")
    for vert in vertices:
        vertfile.write("{},{},{}\n".format(vert[0], vert[1], vert[2]))
    vertfile.close()

    args = [
        multivalue_bin,
        filename_base + ".csv",
        filename_base + ".dx",
        filename_base + "_out.csv",
    ]
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=directory)
    stdout, stderr = p2.communicate()

    # Read the charge file
    chargefile = open(tmp_file_base + "_out.csv")
    charges = numpy.array([0.0] * len(vertices))
    for ix, line in enumerate(chargefile.readlines()):
        charges[ix] = float(line.split(",")[3])

    remove_fn = os.path.join(directory, filename_base)
    os.remove(remove_fn)
    os.remove(remove_fn+'.csv')
    os.remove(remove_fn+'.dx')
    os.remove(remove_fn+'.in')
    os.remove(remove_fn+'-input.p')
    os.remove(remove_fn+'_out.csv')

    return charges

def computeAPBS(vertices, path_pqr, model):
    """
        Calls APBS, pdb2pqr, and multivalue and returns the charges per vertex
    """

    directory = str(PATH_BASE.joinpath("tempfiles"))
    path_apbs = rewrite_apbs_input(path_pqr, model)
    args = [apbs_bin, path_apbs]

    print(args)
    p2 = Popen(args, stdout=PIPE, stderr=PIPE,cwd=directory)
    stdout, stderr = p2.communicate()

    vertpath = Path(directory).joinpath("test.csv")
    with open(vertpath, "w") as vertfile:
        for vert in vertices:
            vertfile.write("{},{},{}\n".format(vert[0], vert[1], vert[2]))

    path_apbs_out = Path(directory).joinpath(f"{path_pqr.name}.dx")
    path_multi = Path(directory).joinpath(f"{path_pqr.name}_out.csv")
    args = [
        multivalue_bin,
        vertpath,
        path_apbs_out,
        path_multi
    ]

    print(args)
    p2 = Popen(args, stdout=PIPE, stderr=PIPE, cwd=directory)
    stdout, stderr = p2.communicate()

    # Read the charge file
    with open(path_multi, 'r') as chargefile:
        charges = numpy.array([0.0] * len(vertices))
        for ix, line in enumerate(chargefile.readlines()):
            charges[ix] = float(line.split(",")[3])

#   remove_fn = os.path.join(directory, filename_base)
    #os.remove(remove_fn)
#   os.remove(remove_fn+'.csv')
#   os.remove(remove_fn+'.pqr.dx')
#   os.remove(remove_fn+'.pdb'+'.in')
#   #os.remove(remove_fn+'-input.p')
#   os.remove(remove_fn+'_out.csv')

    return charges

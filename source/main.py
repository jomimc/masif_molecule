#!/usr/bin/python
import argparse
from pathlib import Path
import os
import shutil
import sys

import numpy as np
import Bio
from Bio.PDB import * 

# Local includes
from default_config.masif_opts import masif_opts
from triangulation.computeMSMS import computeMSMS
from triangulation.fixmesh import fix_mesh
import pymesh
from input_output.extractPDB import extractPDB
from input_output.save_ply import save_ply
from input_output.read_ply import read_ply
from input_output.protonate import protonate
from triangulation.computeHydrophobicity import computeHydrophobicity
from triangulation.computeCharges import computeCharges, assignChargesToNewMesh
from triangulation.computeAPBS import computeAPBS
from triangulation.compute_normal import compute_normal
from sklearn.neighbors import KDTree

if len(sys.argv) <= 1: 
    print("Usage: {config} "+sys.argv[0]+" PDBID_A")
    print("A or AB are the chains to include in this surface.")
    sys.exit(1)



def parse_arguments():
    parser = argparse.ArgumentParser(description='A program to compute molecular surfaces')
    parser.add_argument('path', type=Path, help='Path to molecule PDB / PQR file')
    parser.add_argument('-c', '--chain', default='', help='Choose a single chain to compute')
    parser.add_argument('-m', '--method', type=str, help='Which method to run?')
    parser.add_argument('-v', '--verbose', action="store_false", help='Turn off verbose output')
    parser.add_argument('--msms_density', type=float, default=3.0, help='Density of surface triangulation')
    parser.add_argument('--msms_hdensity', type=float, default=3.0, help='Density of surface triangulation')
    parser.add_argument('--msms_probe', type=float, default=1.5, help='Surface triangulation probe radius')
    parser.add_argument('--redo', action='store_true', help='Overwrite surface comparison results')
    parser.add_argument('--noH', action='store_true', help='Do not protonate PDB file?')
    parser.add_argument('--hbond', action='store_true', help='Calculate hydrogen-bonding potential')
    parser.add_argument('--hphob', action='store_true', help='Calculate Kyte-Doolittle hydrophobicity')
    parser.add_argument('--no_apbs', action='store_false', help='Calculate electrostatic potential')
    parser.add_argument('--run_all', action='store_true', help='Run on all surfaces')
    parser.add_argument('--use_mp', action='store_true', help='Use multiprocessing')

    return parser.parse_args()


def compute_surface(args):

    tmp_dir= Path(masif_opts['tmp_dir'])
    main_path = args.path

    if args.chain != '':
        # Extract a single chain
        path_pdb_chain = tmp_dir.joinpath(f"{main_path.stem}_{args.chain}.{main_path.suffix}")
        extractPDB(main_path, path_pdb_chain, args.chain)
        main_path = path_pdb_chain

    if not args.noH:
        # Add hydrogens
        path_protonated = tmp_dir.joinpath(f"{main_path.stem}_H.{main_path.suffix}")
        protonate(str(main_path), str(path_protonated))
        main_path = path_protonated

    # Create xyzrn file

    # Compute MSMS of surface w/hydrogens, 
    msms_args = ["-density", str(args.msms_density), "-hdensity", str(args.msms_hdensity),
                 "-probe", str(args.msms_probe)]
    vertices, faces, normals, names, areas = computeMSMS(args.path, msms_args)

    # Compute "charged" vertices
    if args.hbond:
        vertex_hbond = computeCharges(out_filename1, vertices, names)

    # For each surface residue, assign the hydrophobicity of its amino acid. 
    if args.hphob:
        vertex_hphobicity = computeHydrophobicity(names)


    # Regularize the mesh
    mesh = pymesh.form_mesh(vertices, faces)
    regular_mesh = fix_mesh(mesh, masif_opts['mesh_res'])

    # Compute the normals
    vertex_normal = compute_normal(regular_mesh.vertices, regular_mesh.faces)

    # Assign hbond values to new mesh
    if args.hbond:
        vertex_hbond = assignChargesToNewMesh(regular_mesh.vertices, vertices,\
            vertex_hbond, masif_opts)

    # Assign hphob values to new mesh
    if args.hphob:
        vertex_hphobicity = assignChargesToNewMesh(regular_mesh.vertices, vertices,\
            vertex_hphobicity, masif_opts)

    # Compute the normals
    if not args.no_apbs:
        vertex_charges = computeAPBS(regular_mesh.vertices, out_filename1+".pdb", out_filename1)

    # Convert to ply and save.
    save_ply(out_filename1+".ply", regular_mesh.vertices,\
                        regular_mesh.faces, normals=vertex_normal, charges=vertex_charges,\
                        normalize_charges=True, hbond=vertex_hbond, hphob=vertex_hphobicity)

    if not os.path.exists(masif_opts['ply_chain_dir']):
        os.makedirs(masif_opts['ply_chain_dir'])
    if not os.path.exists(masif_opts['pdb_chain_dir']):
        os.makedirs(masif_opts['pdb_chain_dir'])
    shutil.copy(out_filename1+'.ply', masif_opts['ply_chain_dir']) 
    shutil.copy(out_filename1+'.pdb', masif_opts['pdb_chain_dir']) 

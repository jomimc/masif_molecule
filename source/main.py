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
    parser.add_argument('-o', '--output_dir', default=masif_opts['ply_dir'], help='Output directory path')
    parser.add_argument('-m', '--method', type=str, help='Which method to run?')
    parser.add_argument('-v', '--verbose', action="store_false", help='Turn off verbose output')
    parser.add_argument('--msms_density', type=float, default=3.0, help='Density of surface triangulation')
    parser.add_argument('--msms_hdensity', type=float, default=3.0, help='Density of surface triangulation')
    parser.add_argument('--msms_probe', type=float, default=1.5, help='Surface triangulation probe radius')
    parser.add_argument('--mesh_res', type=float, default=1.0, help='Surface triangulation probe radius')
    parser.add_argument('--redo', action='store_true', help='Overwrite surface comparison results')
    parser.add_argument('--noH', action='store_true', help='Do not protonate PDB file?')
    parser.add_argument('--hbond', action='store_true', help='Calculate hydrogen-bonding potential')
    parser.add_argument('--hphob', action='store_true', help='Calculate Kyte-Doolittle hydrophobicity')
    parser.add_argument('--no_apbs', action='store_false', help='Calculate electrostatic potential')
    parser.add_argument('--run_all', action='store_true', help='Run on all surfaces')
    parser.add_argument('--use_mp', action='store_true', help='Use multiprocessing')

    return parser.parse_args()


def unpack_names(names, features):
    chain, residx, resname, atomtype = [], [], [], []
    names = ['chain', 'residx', 'resname', 'atomtype']
    for n in names:
        c, ri, rn, a = n.split('_')
        chain.append(c)
        residx.append(int(ri))
        resname.append(rn)
        atomtype.append(a)
    for n, f in zip(names, [chain, residx, resname, atomtype]):
        features[n] = f
    return features


def save_features(path, mesh, features):
    # Set precision of data types
    dtype = {'chain':'U1', 'residx':np.int16, 'resname':'U3', 'atomtype':'U1',
             'hbond':np.float16, 'hphob':np.float16, 'charge':np.float16,
             'mean_curvature':np.float16, 'gaussian_curvature':np.float16}

    # Save each feature to a separate file
    for name, feat in features.items():
        p = path.with_name(f"{path.stem}_{name}")
        np.save(p, feat.astype(dtype[name]))


def compute_surface(args):

    tmp_dir= Path(masif_opts['tmp_dir'])
    main_path = args.path
    path_out = Path(args.output_dir).joinpath(f"{args.path.stem}")
    path_ply = path_out.joinpath(f"{args.path.stem}.ply")
    path_feat = path_out.joinpath(f"{args.path.stem}.npy")
    features = {}

    if args.chain != '':
        # Extract a single chain
        path_pdb_chain = tmp_dir.joinpath(f"{main_path.stem}_{args.chain}.{main_path.suffix}")
        extractPDB(main_path, path_pdb_chain, args.chain)
        main_path = path_pdb_chain

    if args.noH:
        # Remove hydrogens
        path_deprotonated = tmp_dir.joinpath(f"{main_path.stem}_noH.{main_path.suffix}")
        deprotonate(str(main_path), str(path_deprotonated))
        main_path = path_deprotonated
    else:
        # Add hydrogens
        path_protonated = tmp_dir.joinpath(f"{main_path.stem}_H.{main_path.suffix}")
        protonate(str(main_path), str(path_protonated))
        main_path = path_protonated

    # Compute MSMS of surface
    msms_args = ["-density", str(args.msms_density), "-hdensity", str(args.msms_hdensity),
                 "-probe", str(args.msms_probe)]
    vertices, faces, normals, names, areas = computeMSMS(str(main_path), msms_args)

    # Compute "charged" vertices
    if args.hbond:
        vertex_hbond = computeCharges(main_path, vertices, names)

    # For each surface residue, assign the hydrophobicity of its amino acid. 
    if args.hphob:
        vertex_hphobicity = computeHydrophobicity(names)


    # Regularize the mesh
    mesh = fix_mesh(pymesh.form_mesh(vertices, faces), args.mesh_res)

    # Compute the normals
    vertex_normal = compute_normal(mesh.vertices, mesh.faces)

    # Find nearest neighbors between old and new mesh
    kdt = KDTree(vertices)
    dists, result = kdt.query(testset, k=4)

    # Assign names (vertex atom/residue info) to new mesh
    names = names[result[:,0]]
    features = unpack_names(names, features)

    # Assign hbond values to new mesh
    if args.hbond:
        vertex_hbond = assignChargesToNewMesh(mesh.vertices, vertices,\
            vertex_hbond, masif_opts, dists=dists, result=result)
        features['hbond'] = vertex_hbond

    # Assign hphob values to new mesh
    if args.hphob:
        vertex_hphobicity = assignChargesToNewMesh(mesh.vertices, vertices,\
            vertex_hphobicity, masif_opts, dists=dists, result=result)
        features['hphob'] = vertex_hphob

    # Compute the normals
    if not args.no_apbs:
        vertex_charges = computeAPBS(mesh.vertices, out_filename1+".pdb", out_filename1)
        features['charge'] = vertex_charges / 10


    # Add curvature to mesh
    mesh.add_attribute("vertex_mean_curvature")
    mesh.add_attribute("vertex_gaussian_curvature")
    features['mean_curvature'] = mesh.get_attribute("vertex_mean_curvature")
    features['gaussian_curvature'] = mesh.get_attribute("vertex_gaussian_curvature")

    # Save mesh and features
    path_ply.parent.mkdir(exist_ok=True)
    save_ply(str(path_ply), mesh)
    save_features(path_feat, features)

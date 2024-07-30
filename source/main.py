#!/usr/bin/python
import argparse
from pathlib import Path

import numpy as np
from Bio.Data.IUPACData import protein_letters_3to1
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
from masif_modules.read_data_from_surface import get_patches_from_surface
from triangulation.computeHydrophobicity import computeHydrophobicity
from triangulation.computeCharges import computeCharges, assignChargesToNewMesh
from triangulation.computeAPBS import computeAPBS
from triangulation.compute_normal import compute_normal
from sklearn.neighbors import KDTree


def parse_arguments():
    parser = argparse.ArgumentParser(description='A program to compute molecular surfaces')
    parser.add_argument('path', type=Path, help='Path to molecule PDB / PQR file')
    parser.add_argument('-c', '--chain', default='', help='Choose a single chain to compute')
    parser.add_argument('-o', '--output_dir', default=masif_opts['ply_dir'], help='Output directory path')
    parser.add_argument('--msms_density', type=float, default=3.0, help='Density of surface triangulation')
    parser.add_argument('--msms_hdensity', type=float, default=3.0, help='Density of surface triangulation')
    parser.add_argument('--msms_probe', type=float, default=1.5, help='Surface triangulation probe radius')
    parser.add_argument('--mesh_res', type=float, default=1.0, help='Surface triangulation probe radius')
    parser.add_argument('--patch_max_dist', type=float, default=9.0, help='Geodesic patch radius')
    parser.add_argument('--patch_max_size', type=int, default=100, help='Maximum number of vertices in patch')
    parser.add_argument('--redo', action='store_true', help='Overwrite surface comparison results')
    parser.add_argument('--noH', action='store_true', help='Do not protonate PDB file?')
    parser.add_argument('--hbond', action='store_true', help='Calculate hydrogen-bonding potential')
    parser.add_argument('--hphob', action='store_true', help='Calculate Kyte-Doolittle hydrophobicity')
    parser.add_argument('--no_apbs', action='store_true', help='Calculate electrostatic potential')
    parser.add_argument('--patches', action='store_true', help='Decompose surface into patches')

    return parser.parse_args()


def unpack_names(names, features):
    chain, residx, resname, atomtype = [], [], [], []
    featnames = ['chain', 'residx', 'resname', 'atomtype']
    for n in names:
        c, ri, rn, a = n.split('_')[:4]
        chain.append(c)
        residx.append(int(ri))
        resname.append(protein_letters_3to1.get(rn, 'X'))
        atomtype.append(a)
    for n, f in zip(featnames, [chain, residx, resname, atomtype]):
        features[n] = f
    return features


def save_features(path, features):
    # Set precision of data types
    dtype = {'chain':'S1', 'residx':np.int16, 'resname':'S1', 'atomtype':'S1',
             'hbond':np.float16, 'hphob':np.float16, 'charge':np.float16,
             'mean_curvature':np.float16, 'gaussian_curvature':np.float16}

    # Save each feature to a separate file
    for name, feat in features.items():
        p = path.with_name(f"{path.stem}_{name}")
        np.save(p, np.array(feat, dtype[name]))


def get_patches(path, mesh, vertex_normals, features, patch_params):
    # Get patches
    input_feat, rho, theta, mask, neigh_idx = get_patches_from_surface(mesh, vertex_normals, features, patch_params)

    # Save patches
    np.save(path.with_name(f"{path.stem}_rho_wrt_center"), rho)
    np.save(path.with_name(f"{path.stem}_theta_wrt_center"), theta)
    np.save(path.with_name(f"{path.stem}_input_feat"), input_feat)
    np.save(path.with_name(f"{path.stem}_mask"), mask)
    np.save(path.with_name(f"{path.stem}_list_indices"), neigh_idx)

    # Save x, y, z
    np.save(path.with_name(f"{path.stem}_X.npy"), mesh.vertices[:,0])
    np.save(path.with_name(f"{path.stem}_Y.npy"), mesh.vertices[:,1])
    np.save(path.with_name(f"{path.stem}_Z.npy"), mesh.vertices[:,2])


def compute_surface(args):

    # Path to directory for temporary files
    tmp_dir= Path(masif_opts['tmp_dir'])

    # Path to input file
    main_path = args.path

    # Path to output folder
    path_out = Path(args.output_dir).joinpath(f"{args.path.stem}")
    path_out.mkdir(parents=True, exist_ok=True)

    # Path for output mesh
    path_ply = path_out.joinpath(f"{args.path.stem}.ply")

    # Path template for features (and patches)
    path_feat = path_out.joinpath(f"{args.path.stem}.npy")

    # Container for features
    features = {}

    # Rewrite PDB file; extract a single chain if needed
    ctxt = '' if args.chain == '' else f'_{args.chain}'
    path_pdb_chain = tmp_dir.joinpath(f"{main_path.stem}{ctxt}{main_path.suffix}")
    extractPDB(main_path, path_pdb_chain, args.chain)
    main_path = path_pdb_chain

    if args.noH:
        # Remove hydrogens
        path_deprotonated = tmp_dir.joinpath(f"{main_path.stem}{main_path.suffix}")
        deprotonate(str(main_path), str(path_deprotonated))
        main_path = path_deprotonated
    else:
        # Add hydrogens
        path_protonated = tmp_dir.joinpath(f"{main_path.stem}{main_path.suffix}")
        protonate(str(main_path), str(path_protonated))
        main_path = path_protonated

    # Compute MSMS of surface
    msms_args = ["-density", str(args.msms_density), "-hdensity", str(args.msms_hdensity),
                 "-probe", str(args.msms_probe)]
    vertices, faces, normals, names, areas = computeMSMS(main_path, msms_args)

    # Compute "charged" vertices
    if args.hbond:
        vertex_hbond = computeCharges(main_path, vertices, names)

    # For each surface residue, assign the hydrophobicity of its amino acid. 
    if args.hphob:
        vertex_hphob = computeHydrophobicity(names)


    # Regularize the mesh
    mesh = fix_mesh(pymesh.form_mesh(vertices, faces), args.mesh_res)

    # Compute the normals
    vertex_normals = compute_normal(mesh.vertices, mesh.faces)

    # Find nearest neighbors between old and new mesh
    kdt = KDTree(vertices)
    dists, result = kdt.query(mesh.vertices, k=4)

    # Assign names (vertex atom/residue info) to new mesh
    names = np.array(names)[result[:,0]]
    features = unpack_names(names, features)

    # Assign hbond values to new mesh
    if args.hbond:
        vertex_hbond = assignChargesToNewMesh(mesh.vertices, vertices,\
            vertex_hbond, masif_opts, dists=dists, result=result)
        features['hbond'] = vertex_hbond

    # Assign hphob values to new mesh
    if args.hphob:
        vertex_hphob = assignChargesToNewMesh(mesh.vertices, vertices,\
            vertex_hphob, masif_opts, dists=dists, result=result)
        features['hphob'] = vertex_hphob

    # Compute the surface charge
    if not args.no_apbs:
        vertex_charges = computeAPBS(mesh.vertices, main_path, tmp_dir)
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

    # Decompose surface into patches
    if args.patches:
        patch_params = {'max_distance':args.patch_max_dist, 'max_shape_size':args.patch_max_size}

        # Add vertices to mesh
        mesh.add_attribute("vertex_nx")
        mesh.add_attribute("vertex_ny")
        mesh.add_attribute("vertex_nz")
        mesh.set_attribute("vertex_nx", vertex_normals[:,0])
        mesh.set_attribute("vertex_ny", vertex_normals[:,1])
        mesh.set_attribute("vertex_nz", vertex_normals[:,2])

        get_patches(path_feat, mesh, vertex_normals, features, patch_params)


if __name__ == "__main__":
    compute_surface(parse_arguments())


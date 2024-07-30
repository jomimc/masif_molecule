from Bio.PDB import *
from default_config.chemistry import radii, polarHydrogens

"""
xyzrn.py: Read a pdb file and output it is in xyzrn for use in MSMS
Pablo Gainza - LPDI STI EPFL 2019
This file is part of MaSIF.
Released under an Apache License 2.0
"""

def output_pdb_as_xyzrn(path_input, path_xyzrn):
    """
        pdbfilename: input pdb / pqr filename
        xyzrnfilename: output in xyzrn format.
    """
    parser = PDBParser()
    struct = parser.get_structure('', path_input)
    outfile = open(path_xyzrn, "w")
    for atom in struct.get_atoms():
        name = atom.get_name()
        residue = atom.get_parent()
        # Ignore hetatms.
        if residue.get_id()[0] != " ":
            continue
        resname = residue.get_resname()
        residx = residue.get_id()[1]
        chain = residue.get_parent().get_id()
        atomtype = name[0]

        if path_input.suffix == '.pqr':
            R = str(atom.get_bfactor())
        elif path_input.suffix == '.pdb':
            if atomtype not in radii:
                continue
            R = radii[atomtype]

        coords = "{:.06f} {:.06f} {:.06f}".format(*atom.get_coord())
        full_id = f"{chain}_{residx:d}_{resname}_{atomtype}_{name}"

        outfile.write(coords + " " + R + " 1 " + full_id + "\n")



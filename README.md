## MaSIF-Molecule :: A reimplementation of the MaSIF pre-processing pipeline for use on arbitrary molecules.

Code used in [Statistical Survey of Chemical and Geometric Patterns on Protein Surfaces as a Blueprint for Protein-mimicking Nanoparticles](https://www.biorxiv.org/content/10.1101/2024.07.18.604221v2).

See [MaSIF](https://github.com/LPDI-EPFL/masif) for the source repository.

## Description

The MaSIF pipeline uses multiple software packages to extract a solvent-accessible
surface, with chemical features such as electrostatic potential or hydrophobicity.
This pipeline consists of several stages:
 * Convert a PDB file to an xyzr file (atomic coordinates and radii).
 * Compute the surface from an xyzr file using [MSMS](http://mgltools.scripps.edu/packages/MSMS/).
 * Regularize the surface using [PyMesh](https://github.com/PyMesh/PyMesh).
 * Convert the PDB file to a PQR (Q:charge, R:atomic radius) file using [PDB2PQR](http://www.poissonboltzmann.org/).
 * Compute the electrostatic potential on the surface using multivalue, and [APBS](http://www.poissonboltzmann.org/).
 * Compute hydrophobicity - mapping amino acid identity to a hydrophobicity score.
 * Compute hydrogen-bonding potential - requires knowledge of amino acid identity.
 * Compute surface curvature

This pipeline cannot work with non-protein molecules due to the key step of assigning
charges when converting files from PDB to PQR. The MaSIF-Molecule allows starting from PQR files,
which are interconvertible with most molecule dynamics trajectory and topology files using
[MDAnalysis](https://docs.mdanalysis.org/2.0.0/documentation_pages/coordinates/PQR.html).
It is impossible to compute hydrophobicity and hydrogen-bonding potential for non-protein
molecules, since these methods only work for amino-acids.

## Software prerequisites 
MaSIF relies on external software/libraries to handle protein databank files and surface files, 
to compute chemical/geometric features and coordinates, and to perform neural network calculations. 
The following is the list of required libraries and programs, as well as the version on which it was tested (in parenthesis).
* [Python](https://www.python.org/) (3.6)
* [reduce](http://kinemage.biochem.duke.edu/software/reduce.php) (3.23). To add protons to proteins. 
* [MSMS](http://mgltools.scripps.edu/packages/MSMS/) (2.6.1). To compute the surface of proteins. 
* [BioPython](https://github.com/biopython/biopython) (1.66) . To parse PDB files. 
* [PyMesh](https://github.com/PyMesh/PyMesh) (0.1.14). To handle ply surface files, attributes, and to regularize meshes.
* PDB2PQR (2.1.1), multivalue, and [APBS](http://www.poissonboltzmann.org/) (1.5). These programs are necessary to compute electrostatics charges.
 
Alternatively you can use the Docker version, which is the easiest to install (See [Docker container](#Docker-container)).
This is **highly recommended** since installing PyMesh is about as painful as unesthetized dental extraction. 

## Installation 
After preinstalling dependencies, add the following environment variables to your path, changing the appropriate directories:

```
export APBS_BIN=/path/to/apbs/APBS-1.5-linux64/bin/apbs
export MULTIVALUE_BIN=/path/to/apbs/APBS-1.5-linux64/share/apbs/tools/bin/multivalue
export PDB2PQR_BIN=/path/to/apbs/apbs/pdb2pqr-linux-bin64-2.1.1/pdb2pqr
export PATH=$PATH:/path/to/reduce/
export REDUCE_HET_DICT=/path/to/reduce/reduce_wwPDB_het_dict.txt
export PYMESH_PATH=/path/to/PyMesh
export MSMS_BIN=/path/to/msms/msms
export PDB2XYZRN=/path/to/msms/pdb_to_xyzrn
```

Clone masif to a local directory

```
git clone https://github.com/jomimc/masif_molecule
cd masif_molecule/
```

Since MaSIF is written in Python, no compilation is required.

## Docker container

The easiest way to test MaSIF is through a Docker container. Please see our tutorial on reproducing the paper results here:

[Docker container](docker_tutorial.md)


## License

MaSIF-molecule is released under an [Apache v2.0 license](LICENSE).

## Reference
If you use this code, please use the bibtex entry in [citation.bib](citation.bib)

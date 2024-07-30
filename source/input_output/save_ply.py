import pymesh
import numpy
"""
read_ply.py: Save a ply file to disk using pymesh and load the attributes used by MaSIF. 
Pablo Gainza - LPDI STI EPFL 2019
Released under an Apache License 2.0
"""


def save_ply(path, mesh):
    """ Save vertices, mesh in ply format.
    """
    pymesh.save_mesh(path, mesh, use_float=True, ascii=False)


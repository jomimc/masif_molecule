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
#   if normals is not None:
#       n1 = normals[:, 0]
#       n2 = normals[:, 1]
#       n3 = normals[:, 2]
#       mesh.add_attribute("vertex_nx")
#       mesh.set_attribute("vertex_nx", n1)
#       mesh.add_attribute("vertex_ny")
#       mesh.set_attribute("vertex_ny", n2)
#       mesh.add_attribute("vertex_nz")
#       mesh.set_attribute("vertex_nz", n3)

    pymesh.save_mesh(path, mesh, use_float=True, ascii=False)


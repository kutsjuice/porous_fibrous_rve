# -*- coding: utf-8 -*-
"""
Created on Sun Feb 12 13:52:20 2023

@author: Mikhail Kuts
"""

import gmsh
import math
import sys
import numpy as np

gmsh.initialize()

gmsh.model.add("fibers_test")

btm_bnd = -2.5
top_bnd = 2.5
lft_bnd = 2.5
rgh_bnd = -2.5


cyl1 = gmsh.model.occ.addCylinder(0, btm_bnd, 0, 0, top_bnd-btm_bnd, 0, 1)

cyl2 = gmsh.model.occ.addCylinder(lft_bnd, 0, 0, rgh_bnd - lft_bnd, 0, 0, 1)
R = 1.0
Rs = 0.7

union = gmsh.model.occ.fuse([(3,cyl1)], [(3, cyl2)], 3)
# gmsh.model.occ.addCylinder(-2 * R, 0, 0, 4 * R, 0, 0, Rs, 4)
# gmsh.model.occ.addCylinder(0, -2 * R, 0, 0, 4 * R, 0, Rs, 5)
# gmsh.model.occ.addCylinder(0, 0, -2 * R, 0, 0, 4 * R, Rs, 6)
# gmsh.model.occ.fuse([(3, 4), (3, 5)], [(3, 6)], 7)

# gmsh.model.occ.fuse([(3, 1)], [(3, 2)], 7)




gmsh.model.occ.synchronize()

# gmsh.model.addPhysicalGroup(3, [union])


# gmsh.model.mesh.setSize(gmsh.model.getEntities(0), 1.0)

# boundary = gmsh.model.getBoundary([(3,3)], False, False, True)

# gmsh.model.mesh.setSize(gmsh.model.getBoundary([], False, False, True),
#                         lcar3)

# gmsh.model.mesh.setSize(boundary, size=0.2)

gmsh.model.mesh.generate(3)

gmsh.write("test_fiber.vtk")

gmsh.finalize()


# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 12:25:09 2023

@author: Mikhail Kuts
"""
import gmsh
import numpy as np
from numpy.linalg import norm
import SpatialFiber as sf


diam = 0.05


gmsh.initialize()

gmsh.model.add("fibers_test")

fib_a = sf.Str8Fiber([-2.742606148578041192e-01, 
                      -4.394095575269409970e-01,
                      -1.024999999999999911e+00], 
                     [-1.514672052753974207e-01,
                      -2.728893913042389485e-01,
                      1.024999999999999911e+00], diam)

fib_b = sf.Str8Fiber([7.604076057532929145e-01,
                      6.129342736084144594e-01,
                      -1.024999999999999911e+00], 
                     [6.883149987237053846e-01,
                      5.843653078673829437e-01,
                      1.024999999999999911e+00], diam)

print(sf.sk_min_dist(fib_a, fib_b, ))

f_num_a = fib_a.body_ind
f_num_b = fib_b.body_ind

gmsh.model.occ.synchronize()


# union = gmsh.model.occ.fuse([(3,f_num_a)], [(3, f_num_b)], 3)

# gmsh.model.occ.synchronize()

# bndr = gmsh.model.getBoundary([(3,3)], oriented = False)
# edgs = gmsh.model.getBoundary(bndr, False, False)
# eps = 0.1
# box = [2.5 - diam/2 - eps,
#        2.5 - diam/2 - eps,
#        -diam/2-eps,
#        2.5 + diam/2 + eps,
#        2.5 + diam/2 + eps,
#        1 + diam/2 + eps]

# box_edges = gmsh.model.occ.getEntitiesInBoundingBox(*box, dim = 1)

# gmsh.model.occ.fillet([3], [abs(i[1]) for i in box_edges], [0.05])

gmsh.model.occ.synchronize()


gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 10)
gmsh.model.mesh.generate(3)

gmsh.write("check_intersection.vtk")

gmsh.finalize()
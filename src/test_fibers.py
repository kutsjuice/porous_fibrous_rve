# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 12:25:09 2023

@author: Mikhail Kuts
"""
import gmsh
import numpy as np
from numpy.linalg import norm
import SpatialFiber as sf


diam = 0.2


gmsh.initialize()

gmsh.model.add("fibers_test")

fib_a = sf.Str8Fiber([2.5,0,0], [2.5,5,0], diam)
fib_b = sf.Str8Fiber([0,2.5,0.1], [5,2.5,0.1], diam)

f_num_a = fib_a.body_ind
f_num_b = fib_b.body_ind

gmsh.model.occ.synchronize()


union = gmsh.model.occ.fuse([(3,f_num_a)], [(3, f_num_b)], 3)

gmsh.model.occ.synchronize()

bndr = gmsh.model.getBoundary([(3,3)], oriented = False)
edgs = gmsh.model.getBoundary(bndr, False, False)
eps = 0.1
box = [2.5 - diam/2 - eps,
       2.5 - diam/2 - eps,
       -diam/2-eps,
       2.5 + diam/2 + eps,
       2.5 + diam/2 + eps,
       1 + diam/2 + eps]

box_edges = gmsh.model.occ.getEntitiesInBoundingBox(*box, dim = 1)

gmsh.model.occ.fillet([3], [abs(i[1]) for i in box_edges], [0.05])

gmsh.model.occ.synchronize()


fib_c = sf.Str8Fiber([12.5, 0, 0], [12.5, 5, 0], diam)
fib_d = sf.Str8Fiber([15.0, 0, 0.1], [10.0, 5, 0.1], diam)

f_num_c = fib_c.body_ind
f_num_d = fib_d.body_ind




gmsh.model.occ.synchronize()

old_bnds = gmsh.model.getBoundary([(3,f_num_c),(3, f_num_d)], oriented = False)
old_edgs = gmsh.model.getBoundary(old_bnds, False, False)

union2 = gmsh.model.occ.fuse([(3,f_num_c)], [(3, f_num_d)], 6)

gmsh.model.occ.synchronize()

new_bnds = gmsh.model.getBoundary([(3,6)], oriented = False)
new_edgs = gmsh.model.getBoundary(new_bnds, False, False)

diff = set(new_edgs).difference(set(old_edgs))

dist, t0, t1 = sf.sk_min_dist(fib_c, fib_d)
pc = fib_c.get_point(t0)
pd = fib_d.get_point(t1)

mid_point = (pc + pd) / 2

# box = sf.calc_box(mid_point, diam)
# box_edges = gmsh.model.occ.getEntitiesInBoundingBox(*box, dim = 1)

gmsh.model.occ.fillet([6], [abs(i[1]) for i in new_edgs], [0.05])

gmsh.model.occ.synchronize()
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 10)
gmsh.model.mesh.generate(3)

gmsh.write("intersected_fibers.vtk")

gmsh.finalize()
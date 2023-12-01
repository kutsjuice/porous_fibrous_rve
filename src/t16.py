# -*- coding: utf-8 -*-
"""
Created on Sat Feb 11 13:55:45 2023

@author: MKuts
"""

import gmsh
import math
import sys
import numpy as np

gmsh.initialize()

gmsh.model.add("t16")

gmsh.logger.start()

box1 = gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
box2 = gmsh.model.occ.addBox(0, 0, 0, 0.5, 0.5, 0.5)

result_body = gmsh.model.occ.cut([(3, box1)], [(3, box2)],3)

center = np.array([0,0.75,0])

radius = 0.09

holes = []

for i in range(1,6):
    center[[0,2]] += 0.166
    t = gmsh.model.occ.addSphere(*center, radius,)
    holes.append((3,t))

ov, ovv = gmsh.model.occ.fragment([(3,3)], holes)

print("fragment produced volumes:")
for e in  ov:
    print(e)


print("before?after fragment relation:")
for e in zip([(3,3)] + holes, ovv):
    print("parent" + str(e[0]) + " -> child " + str(e[1]))
    


gmsh.model.occ.synchronize()


for i in range(1,6):
     gmsh.model.addPhysicalGroup(3, [3+i], i)

gmsh.model.addPhysicalGroup(3, [ov[-1][1],10])



lcar1 = 0.1
lcar2 = 0.0005
lcar3 = 0.055

gmsh.model.mesh.setSize(gmsh.model.getEntities(0), lcar1)


gmsh.model.mesh.setSize(gmsh.model.getBoundary(holes, False, False, True), lcar3)

eps = 1e-3
ov = gmsh.model.getEntitiesInBoundingBox(0.5 - eps, 0.5 - eps, 0.5 - eps,
                                         0.5 + eps, 0.5 + eps, 0.5 + eps, 0)

gmsh.model.mesh.setSize(ov, lcar2)
gmsh.model.mesh.generate(3)

gmsh.write("t16.vtk")

gmsh.finalize()

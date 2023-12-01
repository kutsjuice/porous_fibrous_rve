# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 21:43:49 2023

@author: MKuts
"""

import gmsh
import numpy as np
from numpy.linalg import norm

from scipy.stats import halfnorm

import SpatialFiber as sf





gmsh.initialize()

gmsh.model.add("fibers_test")

seed = 100

RG = np.random.default_rng(seed)

width = 2
height = 2
depth = 2
diam = 0.15

x_min = -width/2
x_max = width/2

y_min = -height/2
y_max = height/2

z_min = -depth/2
z_max = depth/2

fibers_num = 50

xy = RG.uniform(-1, 1, [fibers_num, 2]) @ np.array([[width/2,0],[0, height/2]])
alpha = RG.uniform(0, 2*np.pi, [fibers_num,1])
beta = np.pi/2 - (halfnorm.rvs(size = [fibers_num, 1]) * 25 * np.pi / 180 + 10 * np.pi / 180) 

fibers = []
fibers_keys = {}

for i in range(fibers_num):
    fibers.append(sf.generate_fiber_by_angles(xy[i,:], alpha[i,0], beta[i,0], [z_min, z_max], diam))
    fibers_keys[i+1] = [fibers[-1]]
k=0
latest_index = len(fibers)
gmsh.model.occ.synchronize()
for i in range(fibers_num):
    
    for j in range(i+1,fibers_num):
        dist, t_a, t_b = sf.sk_min_dist(fibers[i], fibers[j])
        print(k)
        k+=1
        if dist < 0.99*diam:
            gmsh.model.occ.synchronize()
            print("intersected (",i,",",j,")",)
            latest_index += 1
            f1a = fibers[i].start_point
            f1b = fibers[i].end_point
            f2a = fibers[j].start_point
            f2b = fibers[j].end_point
            # print("(",fibers[i].body_ind,",",fibers[j].body_ind,")")
            if fibers[i].body_ind != fibers[j].body_ind:
                buf_ind = gmsh.model.occ.fuse([(3,fibers[i].body_ind)], [(3, fibers[j].body_ind)], latest_index)
                # print("new index = ", latest_index)
                fibers_keys[latest_index] = fibers_keys[fibers[i].body_ind] + fibers_keys[fibers[j].body_ind]
                # print(fibers_keys[latest_index])
                fibers_keys[fibers[i].body_ind] = []
                fibers_keys[fibers[j].body_ind] = []
                for fiber in fibers_keys[latest_index]:
                    fiber.body_ind = latest_index


volumes = []

for key in fibers_keys.keys():
    if fibers_keys[key]:
        volumes.append((3,key))
# sphere = gmsh.model.occ.addSphere(0, 0, 0, 1.5)
            
cutters = []

cutters.append((3, gmsh.model.occ.addBox(4*x_min, 4*y_min, 4*z_min, 3*width/2, 3*height, 3*depth)))
for i in range(3):
    cutters.append(gmsh.model.occ.copy([cutters[-1]])[0])
    rot1 = gmsh.model.occ.rotate([cutters[-1]], 0, 0, 0, 0, 0, 1, np.pi/2)
cutters.append(gmsh.model.occ.copy([cutters[-1]])[0])
gmsh.model.occ.rotate([cutters[-1]], 0, 0, 0, 1, 0, 0, np.pi/2)
cutters.append(gmsh.model.occ.copy([cutters[-1]])[0])
gmsh.model.occ.rotate([cutters[-1]], 0, 0, 0, 1, 0, 0, np.pi)


gmsh.model.occ.cut(volumes, cutters)
gmsh.model.occ.synchronize()
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 5)
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.2)
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.1)
gmsh.option.setNumber("Mesh.MaxNumThreads3D", 4)
gmsh.option.setNumber("Mesh.ElementOrder", 2)
gmsh.option.setNumber("Mesh.Algorithm3D", 9)
gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 2)
# gmsh.option.setNumber("Mesh.RandomFactor", 20)
gmsh.model.mesh.generate(3)

gmsh.write("test_rve.vtk")

gmsh.finalize()
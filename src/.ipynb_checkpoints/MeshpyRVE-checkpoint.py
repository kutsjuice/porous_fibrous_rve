from meshpy.tet import MeshInfo, build
from math import sin, cos, pi, sqrt
import numpy as np
import pyvtk
import random
from Cylinder import cylinder_maker
from merge_meshes import merge_meshes


def fibrous_system(rve_properties, file_name):
    #num_meshes, radius, height, lower_bound, upper_bound, tolerance):
    num_meshes = rve_properties["num_meshes"]
    radius = rve_properties["radius"]
    height = rve_properties["height"]
    width = rve_properties["width"]
    thickness = rve_properties["thickness"]
    center = rve_properties["center"]
    tolerance = rve_properties["tolerance"]
    
    left_bound = center[0] - width/2
    right_bound = center[0] + width/2
    
    lower_bound = center[1] - thickness/2
    upper_bound = center[1] + thickness/2
    
    meshes = [None] * num_meshes
    centers = np.zeros((num_meshes, 3))
    distance_array = np.zeros((1, num_meshes - 1))

    # once the coordinates for the center have been generated, make vector that stores the distance between the new
    # center and all the other centers, and if any entry in the distance vector is less than the radii of both
    # cylinders (plus some tolerance), then a new center will be generated until they no longer intersect

    for i in range(num_meshes):

        # initializing center coordinates, which we will assign to a row in 'centers' array if no intersections occur
        temp_center = [random.uniform(left_bound, right_bound), random.uniform(lower_bound, upper_bound), center[2] - height]

        if i == 0:

            # no intersections can occur in the first iteration, so we just directly assign temp_center
            # to the first row of 'centers', then create the mesh, which is stored in the 'meshes' list

            centers[i, :] = temp_center
            meshes[i] = cylinder_maker(10, height, radius, centers[i, :])

        else:  # after the first iteration, we need to check for intersection

            for j in range(i):
                # we want to compare the distance of temp_center to all other centers, then check for intersections
                distance_array[0, j] = sqrt((temp_center[0] - centers[j, 0]) ** 2 +
                                            (temp_center[1] - centers[j, 1]) ** 2)

            # intersecting is a boolean value that is true if any intersections occur, and false if none occur
            intersecting = np.any((0 < distance_array) & (distance_array < 2 * radius + tolerance))
            # tolerance accounts for non-perfectly cylindrical geometry if few points are used in cylinder_maker

            while intersecting:

                # if intersection occurs with previous temp_center, we want to keep creating new temp_centers until
                # no more intersections occur

                temp_center = [random.uniform(left_bound, right_bound), random.uniform(lower_bound, upper_bound), center[2] - height]

                for j in range(i):

                    distance_array[0, j] = sqrt((temp_center[0] - centers[j, 0]) ** 2 +
                                                (temp_center[1] - centers[j, 1]) ** 2)

                intersecting = np.any((0 < distance_array) & (distance_array < 2 * radius + tolerance))

            # once no intersections are happening, we can now assign temp_center to its appropriate row in 'centers'
            centers[i, :] = temp_center
            meshes[i] = cylinder_maker(10, height, radius, centers[i, :])  # adding the

    # once all the meshes have been stored in the list, we can use merge_meshes to store all of them in a vtk file
    points, cells, bodies = merge_meshes(*meshes)

    vtkelements = pyvtk.VtkData(
        pyvtk.UnstructuredGrid(
            points,
            tetra=cells),
        "Mesh")
    vtkelements.tofile(file_name)
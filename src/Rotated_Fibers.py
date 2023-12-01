import numpy as np
import gmsh

"""
INFO: Two errors occur: "command not done" from the fuse() function and a duplicate facets error can occur after at the
very end of the mesh generation process. Both problems have the same cause, I believe. "intersection could not be 
performed" between fused entity and bounding box can also occur, but is very rare. The script is still totally usable,
but just fails in some ways that are hard to reproduce consistently.
"""


class FiberArray:
    def __init__(self, bounds, center, diameter, length, porosity, tolerance):
        self.bounds = bounds
        self.center = center
        self.diameter = diameter
        self.length = length
        self.porosity = porosity
        self.tolerance = tolerance


# micron scale
width = 50
height = 50
depth = 50

box_dim = np.array([
    [-width, width],
    [-height, height],
    [-depth, depth]
])

diam = 20
leng = 60
target_por = 0.80  # desired void space
tol = 0.01  # tolerance
por = 1
volume_fraction = 0
std = np.pi / 4  # standard deviation of phi angle distribution
mean = 0  # mean of phi angle distribution
itr = 0
move_inside = False

fibers = FiberArray(box_dim, [0, 0, 0], diam, leng, target_por, tol)

gmsh.initialize()


def generate_fiber():  # generates fiber using sphereical coordinates to control orientation and length

    start_coords = -width + (2 * width) * np.random.rand(3)  # xyz coords inside of box
    start_coords = list(start_coords)  # converting np array to list
    # phi = np.random.normal(mean, std)  # gaussian distribution of vertical angle
    phi = np.pi * np.random.rand()  # uniform distribution of vertical angle
    theta = 2 * np.pi * np.random.rand()  # uniform distribution of horizontal angle

    displacement_vec = [None] * 3  # displacement uses spherical coordinates
    displacement_vec[0] = leng * np.sin(phi) * np.cos(theta)  # x-displacement
    displacement_vec[1] = leng * np.sin(phi) * np.sin(theta)  # y-displacement
    displacement_vec[2] = leng * np.cos(phi)  # z-displacement

    gmsh.model.occ.addCylinder(*start_coords, *displacement_vec, diam / 2)

    gmsh.model.occ.synchronize()


gmsh.model.occ.synchronize()
print(gmsh.model.occ.getEntities(3))

coords = []
displacement_list = []
floating = False

# gmsh.option.setNumber("Geometry.AutoCoherence", 2)

while por > target_por:

    """
    This loop sequentially adds fibers until a desired porosity is reached. It's guaranteed that each fiber is
    intersecting at least 1 other fiber, since that's how this fibrous system gains structural integrity. The fibers are 
    fused together into a single body as they are added. Parts of fibers that leave the desired cubic boundaries are
    removed to apply the RVE concept and to enforce boundary conditions. 
    """

    itr += 1
    print("iteration:", itr)

    generate_fiber()

    entities = gmsh.model.occ.getEntities(3)
    if len(entities) == 3:
        print("Unusual output from generate_fiber(), entities are: ", entities)
        gmsh.model.occ.remove([entities[2]])
    elif len(entities) > 3:
        print("Unusual output from generate_fiber(), entities are: ", entities)
        print([entities[2:len(entities)]])
        gmsh.model.occ.remove(entities[2:len(entities)])

    if itr > 1:  # checking intesection, must be after 1st iteration

        print("The entities are:", gmsh.model.getEntities(3))
        intersect_check = gmsh.model.occ.intersect([(3, 1)], [(3, 2)], removeObject=False, removeTool=False)
        gmsh.model.occ.synchronize()
        # Other fibers are fused into single body, so we just need to check if new fiber intersects that one body

        if intersect_check == ([], [[], []]):  # if the intersection is empty
            gmsh.model.occ.remove([(3, 2)], recursive=True)  # removing the fiber if it is floating
            print("Fiber is floating")
            gmsh.model.occ.synchronize()
            floating = True

        # if intersection created two output entities, they both need to be deleted
        elif len(intersect_check[0]) == 1 and intersect_check != ([], [[], []]):
            print("intersected before floating")
            print("pre-intersection-removal ", gmsh.model.occ.getEntities(3))
            intersection_tag = gmsh.model.occ.getEntities(3)[-1][1]  # second element of the last tuple
            gmsh.model.occ.remove([(3, intersection_tag)])  # removing the new intersection entity
            gmsh.model.occ.synchronize()
            print("After intersection-removal ", gmsh.model.occ.getEntities(3))
            gmsh.model.occ.fuse([(3, 1)], [(3, 2)])  # new body tag will be 1, since it is the object
            print("After fusing: ", gmsh.model.occ.getEntities(3))
            gmsh.model.occ.synchronize()
            floating = False

        elif len(intersect_check[0]) >= 2:

            print("intersected after floating")
            print("pre-intersection-removal ", gmsh.model.occ.getEntities(3))
            entities = gmsh.model.getEntities(3)
            intersection_entities = entities[2:len(entities)]
            gmsh.model.occ.remove(intersection_entities)  # removing the new intersection entity
            gmsh.model.occ.synchronize()
            print("After intersection-removal ", gmsh.model.occ.getEntities(3))
            gmsh.model.occ.fuse([(3, 1)], [(3, 2)])  # new body tag will be 1, since it is the object
            gmsh.model.occ.synchronize()
            print("After fusing: ", gmsh.model.occ.getEntities(3))
            gmsh.model.occ.synchronize()
            floating = False

        while floating:
            generate_fiber()

            intersect_check = gmsh.model.occ.intersect([(3, 1)], [(3, 2)], removeObject=False, removeTool=False)
            gmsh.model.occ.synchronize()
            print(intersect_check)
            # Other fibers are fused into single body, so we just need to check if new fiber intersects that one body

            if intersect_check == ([], [[], []]):  # if the intersection is empty
                gmsh.model.occ.remove([(3, 2)], recursive=True)  # removing the fiber if it is floating
                print("Fiber is floating")
                gmsh.model.occ.synchronize()
                floating = True

            elif len(intersect_check[0]) == 1 and intersect_check != ([], [[], []]):
                print("intersected after floating")
                print("pre-intersection-removal ", gmsh.model.occ.getEntities(3))
                intersection_tag = gmsh.model.occ.getEntities(3)[-1][1]  # second element of the last tuple
                gmsh.model.occ.remove([(3, intersection_tag)])  # removing the new intersection entity
                gmsh.model.occ.synchronize()
                print("After intersection-removal ", gmsh.model.occ.getEntities(3))
                gmsh.model.occ.fuse([(3, 1)], [(3, 2)])  # new body tag will be 1, since it is the object
                print("After fusing: ", gmsh.model.occ.getEntities(3))
                gmsh.model.occ.synchronize()
                floating = False

            # if intersection created two or more output entities, they all need to be deleted
            elif len(intersect_check[0]) >= 2:

                print("intersected after floating")
                print("pre-intersection-removal ", gmsh.model.occ.getEntities(3))
                entities = gmsh.model.getEntities(3)
                intersection_entities = entities[2:len(entities)]
                gmsh.model.occ.remove(intersection_entities)  # removing the new intersection entity
                gmsh.model.occ.synchronize()
                print("After intersection-removal ", gmsh.model.occ.getEntities(3))
                gmsh.model.occ.fuse([(3, 1)], [(3, 2)])  # new body tag will be 1, since it is the object
                gmsh.model.occ.synchronize()
                print("After fusing: ", gmsh.model.occ.getEntities(3))
                gmsh.model.occ.synchronize()
                floating = False

    bounding_box = (3, gmsh.model.occ.addBox(-width, -height, -depth, 2 * width, 2 * height, 2 * depth))
    gmsh.model.occ.intersect([(3, 1)], [bounding_box])
    # gmsh.model.occ.fragment([(3, 1)], [bounding_box])
    volume_fraction = gmsh.model.occ.getMass(3, 1)/((box_dim[0, 1] - box_dim[0, 0])**3)
    por = 1 - volume_fraction

    gmsh.model.occ.remove(gmsh.model.occ.getEntities(dim=2))
    gmsh.model.occ.remove(gmsh.model.occ.getEntities(dim=1))
    gmsh.model.occ.remove(gmsh.model.occ.getEntities(dim=0))
    # gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()
    print("after cut: ", gmsh.model.getEntities(3))
    print("Porosity: ", np.round(por, 5))
    # print("Boundary entities: ", gmsh.model.getBoundary([(3, 1)]))


gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 5)
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.4 * diam)
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.1 * diam)  # there needs to be a large number of elements to accurately
gmsh.option.setNumber("Mesh.MaxNumThreads3D", 4)        # mesh the complex surfaces created by intersection
gmsh.option.setNumber("Mesh.MeshSizeFactor", 0.4)
gmsh.option.setNumber("Mesh.Algorithm3D", 9)  # R-tree algorithm, better at handling intersections than Delaunay
gmsh.option.setNumber("Mesh.ElementOrder", 2)  # First order elements

gmsh.model.occ.synchronize()

gmsh.model.mesh.generate(3)
# gmsh.model.mesh.optimize("UntangleMeshGeometry")
# gmsh.model.mesh.recombine()


file_path = 'out_mesh'
gmsh.write(f"{file_path}/por_{int(1000 * np.round(por, 3))}_size_0.4.vtk")  # filename is based on porosity
gmsh.write(f"{file_path}/por_{int(1000 * np.round(por, 3))}_size_0.4.msh")

gmsh.option.setNumber("Mesh.MeshSizeMax", 1. * diam)
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.6 * diam)
gmsh.model.mesh.generate(3)
gmsh.write(f"{file_path}/por_{int(1000 * np.round(por, 3))}_size_1.0.vtk")  # filename is based on porosity
gmsh.write(f"{file_path}/por_{int(1000 * np.round(por, 3))}_size_1.0.msh")

gmsh.option.setNumber("Mesh.MeshSizeMax", 0.1 * diam)
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.05 * diam)
gmsh.model.mesh.generate(3)
gmsh.write(f"{file_path}/por_{int(1000 * np.round(por, 3))}_size_0.1.vtk")  # filename is based on porosity
gmsh.write(f"{file_path}/por_{int(1000 * np.round(por, 3))}_size_0.1.msh")

print("Porosity", np.round(por, 5))
gmsh.finalize()

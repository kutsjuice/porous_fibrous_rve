import gmsh

ORDER = 2
## Setting up gmsh properties
gmsh.initialize()

# Choose if Gmsh output is verbose
gmsh.option.setNumber("General.Terminal", 0)

# Set elements order to the specified one
gmsh.option.setNumber("Mesh.ElementOrder", ORDER)
# Set elements size
# gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 5) # uncomment to use for mesh refinement dependending from its surface curvature
gmsh.option.setNumber("Mesh.MeshSizeMax", 5e-2)
gmsh.option.setNumber("Mesh.MeshSizeMin", 1e-2)

# Set threads number for distrebuted meshing
# gmsh.option.setNumber("Mesh.MaxNumThreads3D", 4)

# Set mesh algorithm (default is Delaunay triangulation)
# see https://gmsh.info/doc/texinfo/gmsh.html#Choosing-the-right-unstructured-algorithm
gmsh.option.setNumber("Mesh.Algorithm3D", 3)

# gmsh.option.setNumber("Mesh.RecombinationAlgorithm",3)
# gmsh.option.setNumber("Mesh.Recombine3DAll",1)

# Set the usage of hexahedron elements 
gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 0)
## Importing RVE geometry
gmsh.open("./out_mesh/por_0.783.msh")

model = gmsh.model()
# model.add("main_domain")

# Synchronize OpenCascade representation with gmsh model
model.occ.synchronize()


# Generate the mesh
# model.mesh.generate(2)
# model.mesh.recombine()
model.mesh.generate(dim=3)

file_path = 'out_mesh'
gmsh.write(f"{file_path}/por_783_01.vtk")  # filename is based on porosity
gmsh.write(f"{file_path}/por_783_01.msh")
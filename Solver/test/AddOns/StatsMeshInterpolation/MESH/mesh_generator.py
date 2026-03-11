import gmsh
import sys
import meshio
import numpy as np

meshFilename = "mesh"

x0 = 0.0
x1 = 1.0
y0 = 0.0
y1 = 2.0
z0 = 0.0
z1 = 3.0

numElem_y = [3]
numElem_z = [2]


geoKernel = gmsh.model.geo

gmsh.initialize()

# Create a model
gmsh.model.add("ExtrudedMesh")

# Create points
lc = 0.5
Px0 = geoKernel.addPoint(x0,y0,z0,lc)
Px1 = geoKernel.addPoint(x1,y0,z0,lc)

# Define line
line1 = geoKernel.addLine(Px0, Px1)
curve1 = geoKernel.addCurveLoop([line1])

# Synchronize
geoKernel.synchronize()

# Mesh the line
gmsh.model.mesh.generate(1)

A = geoKernel.extrude([(1,curve1)], 0, y1-y0, 0, numElements=numElem_y, heights=[1], recombine=True)

bottomSurfTag = [tag for (dim,tag) in A if dim == 2][0]

A = geoKernel.extrude([(2,bottomSurfTag)], 0, 0, z1-z0, numElements=numElem_z, heights=[1], recombine=True)

volTag = [tag for (dim,tag) in A if dim == 3][0]

topSurfTag = -1
leftSurfTag = -1
rightSurfTag = -1
frontSurfTag = -1
backSurfTag = -1

TOLZero = 1e-8

geoKernel.synchronize()

for (dim,tag) in gmsh.model.getEntities(dim=2):
    n = gmsh.model.getNormal(tag, [0.5, 0.5])
    if ( np.linalg.norm(n - np.array([0,0,1], dtype=float)) < TOLZero ):
        topSurfTag = tag
    if ( np.linalg.norm(n - np.array([-1,0,0], dtype=float)) < TOLZero ):
        leftSurfTag = tag
    if ( np.linalg.norm(n - np.array([1,0,0], dtype=float)) < TOLZero ):
        rightSurfTag = tag
    if ( np.linalg.norm(n - np.array([0,-1,0], dtype=float)) < TOLZero ):
        frontSurfTag = tag
    if ( np.linalg.norm(n - np.array([0,1,0], dtype=float)) < TOLZero ):
        backSurfTag = tag

geoKernel.addPhysicalGroup(2, [bottomSurfTag], name="bottom")
geoKernel.addPhysicalGroup(2, [topSurfTag], name="top")
geoKernel.addPhysicalGroup(2, [leftSurfTag], name="left")
geoKernel.addPhysicalGroup(2, [rightSurfTag], name="right")
geoKernel.addPhysicalGroup(2, [frontSurfTag], name="front")
geoKernel.addPhysicalGroup(2, [backSurfTag], name="back")

geoKernel.addPhysicalGroup(3, [volTag], name="vol")

geoKernel.synchronize()
gmsh.model.mesh.generate(3)

# Write mesh to file
gmsh.option.setNumber("Mesh.MshFileVersion", 2)
gmsh.write(meshFilename + ".msh")
print("Mesh written to " + meshFilename)

# Finalize Gmsh
gmsh.finalize()


M = meshio.read(meshFilename + ".msh")
M.write(meshFilename + ".vtk")
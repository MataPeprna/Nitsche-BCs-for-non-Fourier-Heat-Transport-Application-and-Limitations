# This code has been based on a code by Jaroslav Hron from lecture NMMO403 at MFF UK.

import sys, os
import numpy as np
from dolfin import *
import mshr
import math


def py_snap_boundary(mesh, sub_domain):
    boundary = BoundaryMesh(mesh,"exterior")
    dim = mesh.geometry().dim()
    x = boundary.coordinates()
    for i in range(0,boundary.num_vertices()):
        sub_domain.snap(x[i,:])
    ALE.move(mesh,boundary)


# Define domain perematers
center = Point(10, 2.05)
radius = 0.5
refinement_size = 0.1
L = 20
W = 4

# Build domain and mesh
def generate(n=0):
    begin("Generating meshes...")

    # Build mesh
    geometry= mshr.Rectangle(Point(0.0,0.0), Point(L, W)) - mshr.Circle(center, radius, 50) #50
    mesh = mshr.generate_mesh(geometry, 200)

    mesh.init()
    
    class Cylinder(SubDomain):
        def snap(self, x):
            r = sqrt((x[0] - center[0])**2 + (x[1] - center[1])**2)
            if r <= radius:
                x[0] = center[0] + (radius / r)*(x[0] - center[0])
                x[1] = center[1] + (radius / r)*(x[1] - center[1])
        def inside(self, x, on_boundary):
            r = sqrt((x[0] - center[0])**2 + (x[1] - center[1])**2)
            return( (r<=2.0*radius) and on_boundary )
                
    cylinder=Cylinder()

    def mark_bndry(mesh):
        # Construct facet (here boundary edges) markers
        bndry = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
        for f in facets(mesh):
            if f.exterior():
                mp = f.midpoint()
                if near(mp[0], 0.0): bndry[f] = 1                      # Inflow
                elif near(mp[0], L): bndry[f] = 2                      # Outflow
                elif near(mp[1], 0.0) or near(mp[1], W): bndry[f] = 3  # Walls
                elif mp.distance(center) <= radius: bndry[f] = 5       # Cylinder
        return(bndry)

    bndry=mark_bndry(mesh)

    width = refinement_size
    parameters["refinement_algorithm"] = "plaza_with_parent_facets"

    # Refinement around the cylinder (for cells with distance <= refinement-size). Higher levels of refinement are performed on exponentially smaller
    # neighbourhoods of the cylinder to reduce number of cells and speed up the code.
    for k in range(n):
        info("refinement level {}".format(k+1))
        cf = MeshFunction('bool', mesh, mesh.topology().dim(), False)
        treshold = radius + width
        for c in cells(mesh):
            if c.midpoint().distance(center)<treshold  : cf[c]=True  # Local refinement
            #cf[c]=True # Global refinement. Not used.
        width = width/2
            
        mesh=refine(mesh,cf)
        py_snap_boundary(mesh,cylinder)
        mesh.snap_boundary(cylinder, False)
        bndry=mark_bndry(mesh)
        cylinder.mark(bndry, 5)
        
    meshes = [(mesh, bndry)]
    return(meshes)






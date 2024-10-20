# This code has been inspired by a code by Jaroslav Hron for lecture NMMO403 at MFF UK.

import numpy as np
from dolfin import *
import mshr
import math


# Define domain perameters
center = Point(10, 2.05)
radius = 0.5
L = 20
W = 4

# Generate domain
def generate(n=0):
    begin("Generating meshes...")

    # Build mesh
    geometry= mshr.Rectangle(Point(0.0,0.0), Point(L, W)) - mshr.Circle(center, radius, 50) # Circle is approximated as regular 50-gon.
    mesh = mshr.generate_mesh(geometry, 200)                                                

    mesh.init()

    # Construct facet (here boundary edges) markers
    def mark_bndry(mesh):
        bndry = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
        for f in facets(mesh):
            if f.exterior():
                mp = f.midpoint()
                if near(mp[0], 0.0): bndry[f] = 1                          # Inflow
                elif near(mp[0], L): bndry[f] = 2                          # Outflow
                elif near(mp[1], 0.0) or near(mp[1], W): bndry[f] = 3      # Walls
                elif mp.distance(center) <= radius: bndry[f] = 5           # Cylinder
        return(bndry)

    bndry=mark_bndry(mesh)

    # Refinement algorithm (for n =! 0) is non-trivial and written by Jaroslav Hron. In this code, it is not included due to copyright issues.
    # Global refinement can be performed by increasing the number parameter in mshr.generate_mesh above.
    # Increasing the number n times decreases the edge length roughly n-times.
    
    meshes = [(mesh, bndry)]
    return(meshes)


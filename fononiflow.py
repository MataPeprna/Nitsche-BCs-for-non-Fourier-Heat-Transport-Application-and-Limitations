# This code has been inspired by a code by Jaroslav Hron from lecture NMMO403 at MFF UK (not accessible publicly)
# and a similar code accessible at https://fenics-course-fluids.readthedocs.io/en/latest/NavierStokes2D/NavierStokes2D.html

from dolfin import *
import mshr
import numpy as np
import math
import matplotlib.pyplot as plt

import mymeshsimple

# Set how many iterationf of refinement around the cylinder should be done (here 1)
[(mesh, bndry)] = mymeshsimple.generate(0)   # No refinement. In the simplified code omitting the non-trivial code by Jaroslav Hron, this refinement
                                             # is not possible. One can refine the mesh directly in mymeshsimple.py file by changing the
                                             # parameter in mrhr.generate_mesh

# Define finie elements; here Tylor-Hood
Ee = FiniteElement("CG",mesh.ufl_cell(),1) 
Em = VectorElement("CG",mesh.ufl_cell(),2)
Eme=MixedElement([Em,Ee]) 

# Define function spaces
M = FunctionSpace(mesh,Em)
E = FunctionSpace(mesh,Ee)
W = FunctionSpace(mesh,Eme)

# Dirichlet boundary conditions
noslip = Constant((0, 0))
hot = Constant(4.024)  
cold = Constant(3.976)  
bcm_walls = DirichletBC(W.sub(0), noslip, bndry, 3) # Dirichlet BC on walls. Here no-slip.
bcm_cylinder = DirichletBC(W.sub(0), noslip, bndry, 5) # Dirichlet BC on cylinder. Here no-slip. Not used.
bce_left = DirichletBC(W.sub(1), hot, bndry, 1) # BC for energy density e
bce_right = DirichletBC(W.sub(1), cold, bndry, 2) # BC for energy density e


# Parameters for slip BC
Theta = 0.0  # Principal parameter in Navier-slip BC, setting how slippery the boundary is.
             # 0=full slip.
             # Values approaching 1 converge to no slip, but Theta=1 is not well defined because of division by 0. No-slip has to be imposed separately.
             # All values in between work.
gamma = 1    # Another parameter in Navier-slip compensating for pyhsical units. Can be left constant as its change has the same effects as changing Theta. 
beta = 1000  # Parameter beta from stabilisation term in Nitsche method. The higher, the larger the effect of the stabilisation.


# Specify material parameters and time stepping
#en = Constant(1)
tau = Constant(0.00005)   # Dissipation time
taur = Constant(25000/4)  # (Another) Dissipation time 
ce = Constant(8)          # Sound speed

viscosity = tau*ce*ce/5   # (Computed) equivalent to Navier-stokes viscosity 

dt = 20                   # Time step
t_end = 1000              # Final time
theta=Constant(1.0)       # Time stepping (0 = Forward Euler, 0.5 = Crank-Nicolson, 1 = Backward Euler); BE is stable but converges slowly, which does not
                          # matter here as the solution is stationary in time (after an initialisation period).

# Inflow/outflow boundary condition for velocity and initial condition
U=15
w_init = Expression(("(((-0.005)*(1)*cosh(x[1]/2-1)/cosh(1))+0.005*(1))*(x[0]-10)*(x[0]-10)/100","0","4.024 - 0.0024*x[0]"), degree=0) # Specify IC.
m_in = Expression(("((-0.005)*(1)*cosh(x[1]/2-1)/cosh(1))+0.005*(1) ", "0.0"), U=U, degree=2) # This m is a stable solution. 
bcm_in = DirichletBC(W.sub(0), m_in, bndry, 1)   # Specify inflow BC
bcm_out = DirichletBC(W.sub(0), m_in, bndry, 2)  # Specify outflow BC

# Load initial condition from file; not used
#m_init = Function(M, 'saved_m.xml')
#e_init = Function(E, 'saved_e.xml')

# Specify boundary onditions used for the computation
bcs = [bcm_in, bcm_out, bcm_walls, bce_left]

# Define quantities used in the weak formulation of the problem
n = FacetNormal(mesh)                      # Outer normal to facets (in 2D edges of the triangulation)
edgelen = mesh.hmin()                      # Lenght of the longes edge in the triangulation
I = Identity(mesh.geometry().dim())        
ds = Measure("ds", subdomain_data=bndry)   # Define integration over boundary

# Define/name unknown and test function(s)
(m_, e_) = TestFunctions(W)

# Unknown at current time step
w = Function(W)
(m, e) = split(w)

# Already known soution at the previous time step
w0 = Function(W)
(m0, e0) = split(w0)

# Define weak forms

def pressure(e,m, mk, ek):
    return e/3+ce*ce*inner(m, mk)/(4*ek)

def StressT(e,m, mk, ek):
    D = sym(grad(m))
    return -pressure(e,m, mk, ek)*I + 0.2*ce*ce*tau*2*D - (2/15)*ce*ce*tau*div(m)*I

def vn(v,n):
    return inner(v,n)*n

def vt(v,n):
    return v - vn(v,n)

def a(v,u) :
    D = sym(grad(v))
    return 0.2*ce*ce*tau*inner(2*D, grad(u))*dx

def b(q,v) :
    return inner(div(v),q)*dx

def c(v,u) :
    return inner(v/taur,u)*dx

def d(v,u) :
    return inner(div(v),2*ce*ce*tau*div(u)/15)*dx

def f(v,u, e) :
   return 3*ce*ce*(1/(4*e))*inner(v, grad(u)*v)*dx # for our model
    #return 0 # for Guyer-Krumhansl

def N0(m,m_) :
    if Theta == 0.0:
        return 0
    else:
        return (Theta/(gamma*(1.0-Theta))*inner(vt(m,n),vt(m_,n)))*ds(5) 
    
def N1(e,m,m_):
    return (-inner(dot(StressT(e,m,m,e),n),n)*inner(m_,n))*ds(5) 
    
def N2(e_,m_,m,mk, ek):
    return (inner(dot(StressT(e_,m_,mk, ek),n),n)*inner(m,n))*ds(5) 

def Nstab(m,m_):
    return beta*viscosity/edgelen * (inner(m,m_)*ds(5))

# Define pressure
p = e/3 - (ce*ce*dot(m,m))/(4*e)
p0 = e0/3 - (ce*ce*dot(m0,m0))/(4*e0)

# Weak form of the equation without time derivative in current time
F1 = a(m,m_) - b(p,m_) - d(m,m_) + b(e_,ce*ce*m) - f(m,m_,e) + c(m,m_) + N0(m,m_) + N1(e,m,m_) + N2(e_, m_, m, m0, e0) + Nstab(m,m_)

# Weak form of the equation without time derivative in previous time
F0 = a(m0,m_) - b(p0,m_) - d(m0,m_) + b(e_,ce*ce*m0) - f(m0,m_,e) + c(m0,m_) + N0(m0,m_) + N1(e0,m0,m_) + N2(e_, m_, m0, m0, e0) + Nstab(m0,m_)
# For Backward Euler timestepping, this can be omitted. 

# Combine variational forms with time derivative
#  dw/dt + F(w,t) = 0 is approximated as
#  (w-w0)/dt + theta*F(w,t) + (1-theta)*F(w0,t0) = 0

mdot=Constant(1.0/dt)*inner(m-m0,m_)*dx
edot=Constant(1.0/dt)*inner(e-e0,e_)*dx
F = mdot + edot + theta*F1 + (1.0-theta)*F0 

J = derivative(F, w)    # Specifying this derivative accenlerates computation

# Define the problem to be solved (equation, BCs etc.) and the solver to be used
problem=NonlinearVariationalProblem(F,w,bcs,J)
solver=NonlinearVariationalSolver(problem)

assign(w,interpolate(w_init,W))    # Loading initial momentum and energy from an expression within this code.
#assign(w.sub(0),m_init)           # Loading initial momentum from file. Not used.
#assign(w.sub(1),e_init)           # Loading initial energy from file. Not used.

# Specify parameters of the (non-linear) solver
prm = solver.parameters
prm['nonlinear_solver'] = 'newton'
prm['newton_solver']['absolute_tolerance'] = 1E-10
prm['newton_solver']['relative_tolerance'] = 1e-10
prm['newton_solver']['maximum_iterations'] = 20
prm['newton_solver']['linear_solver'] = 'mumps'
prm['newton_solver']['error_on_nonconvergence'] = False

# Create files for storing solution
name="solution_20_200_0_1000" # Named as follows: name _ dt _ n(mesh size) _ refinement level _ beta
out_file=dict()
for i in ['m', 'e'] :
    out_file[i] = XDMFFile(f"results_{name}/{i}.xdmf")
    out_file[i].parameters["flush_output"] = True

(m,e)=w.split(True)
m.rename("m", "momentum")
e.rename("e", "energy")


# Time-stepping
t = 0.0             # Specify initial time.

# Extract individual fields from the vector of all unknowns:
assign(m, w.sub(0))
assign(e, w.sub(1))

# Save to file
out_file['m'].write(m, t)
out_file['e'].write(e, t)

# Initialise field to store postprocessing results
lift=[]
drag=[]
wallcond=[]
incond=[]

# Define functions for postprocessing computations

def through(m,n):
    return inner(m,n)*inner(m,n)

def pres(e,m):
    return e/3+ce*ce*inner(m, m)/(4*e)

def Stress(e,m):
    D = sym(grad(m))
    return -pres(e,m)*I + 0.2*ce*ce*tau*D - (2/15)*ce*ce*tau*div(m)*I

# Compute the solution

while t < t_end:

    # Move current solution to previous slot w0
    w0.assign(w)

    # update time-dependent parameters
    t += dt
    m_in.t=t      # This can be used for time-dependend BCs. Not used.

    # Compute
    begin("Solving ....")
    (_, converged) = solver.solve()
    if not converged:
       break
    end()

    # Extract solutions:
    assign(m, w.sub(0))
    assign(e, w.sub(1))

    # Compute equivalent to Navier-Stokes drag (not used), lift (not used), flux through the cylinder wall (should be zero) and flux through inlet boundary
    #stress = Stress(e,m)
    #force = dot(stress, n)
    #D = -(2.0*force[0]/(1.0*1.0*0.1))*ds(5)
    #L = -(2.0*force[1]/(1.0*1.0*0.1))*ds(5)
    #drag.append((t,assemble(D)))
    #lift.append((t,assemble(L)))
    wallcond.append((t,assemble(through(w.sub(0),n)*ds(5))))
    incond.append((t,assemble(through(w.sub(0),n)*ds(1))))
    
    # Save to file
    out_file['m'].write(m, t)
    out_file['e'].write(e, t)


# Save the computed fluxes into text files
np.savetxt('wallcond_20_200_0_1000.txt',wallcond)
np.savetxt('incond_20_200_0_1000.txt',incond)

# Save the results (to load them as ICs if needed)    
#File('saved_m.xml') << m 
#File('saved_e.xml') << e


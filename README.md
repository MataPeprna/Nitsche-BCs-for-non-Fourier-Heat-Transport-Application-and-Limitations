# Nitsche-Boundary-Conditions-for-non-Fourier-Heat-Transport
Non-Fourier heat conduction is a field of complex models describing heat transport in very specific a challenging settings (i.e. extremely small sizes etc.). This project presents a way to model slip boundary conditions for a complex non-Fourier heat conduction model and implements it numerically. Furthermore, it shows a possible drawback o fhte method and warns about wrong application.

__________________

About this project:

Fourier law states (roughly speaking) that heat travels in the direction, where temperature is decreasing the fastest. This sounds logical and it is usually true. However, in certain complex settings like very small size scale etc., this law is not valid and heat conduction rather resembles flow of a fluid (for instance wind). This sort of heat transport can be described by a rather complex field called phonon hydrodynamics. (Phonons are quasi particles describing for example heat in a crystal and hydrodynamics is a science describing flows od fluids.)

One of the recently developed models of phonon hydrodynamics is particularly similar to the famous Navier-Stokes equations. This projects implements full-slip* (and also Navier-slip) boundary conditions to the model by the so-called Nitsche method. It studies, how to choose a stabilisation parameter to optimise performance, shows that the use of Nitsche method could be dangerous as seeming convergence may occur and suggests further optimisations and check to verify that the computed solutions is indeed correct.

This project contains an in-depth .pdf summary of the results and working codes that give the results from the .pdf (Fononi.py is the main code to run). In order to make the codes work, one needs to download legacy FEniCS 2019.1.0, see [this site](https://fenicsproject.org/download/archive/) for download. Then simply run Python (or Python 3) Fononi.py. For changes in the geometry (lenght of a channel etc.), edit Mymesh_fononi.py. For changes in the level of refinement around the cylindrical obstacle and for changes of material parameters (like relaxation times etc.), edit Fononi.py. This readme then contains the main results without without the technical details in the .pdf document. 

*Full-slip says that there is no flow through the boundary, but the boundary is perfectly slippery and does not slow down the flow in the tangential direction (along itself).

__________________

About the code:

The numerical code is implemented in legacy FEniCS project, version 2019.1.0. It uses backward Euler time stepping (simple first-order method) in space and Taylor-Hood finite element method (FEM) in space (socend order method). Lower order in time is no serious issue as the resulting flow is stable in time. Backward Euler method has the advantage of being stable, reducing shocks caused by the initial conditions.

Nitsche method is a way to enforce boundary conditions weakly, i.e. they are incorporated into the integral forms within the weak FEM formulation. There are two basic types of Nitshce method -- symmetric and nonsymmetric. In general, non-symmetric works without the stabilisation term, while the symmetric needs the stabilisation term with sufficiently large stabilisaiton parameter. Since the non-symmetric variant is more general in the sense above, we use primarily the non-symmetric one.

__________________

About the results:

In short, Nitsche method can be applied to modelling full-slip boundary conditions. The choice of stabilisation parameter does matter. The non-symmetric and symmetric variants seem to give similar results, even though more proper testing of the symmetric Nitsche would be necessary to verify the statement. There are two main problems.

- The method converges extremely slowly (needs extremely fine mesh to converge to the correct result), so it is inapplicable in practical 3D problems.
- A seeming convergence happens, at least in some geometries. For rough meshes, refinement does not change the solution, but the solution is wrong (enforces no-slip boundary ondition instead of full-slip) and only after more refinements, it slowly changes to the correnct full-slip-solution. This seeming convergence can happen on meshes that are seem "reasonably" fine, so there is a danger of actually mistaging the seeming convergence for real convergence.

Finally, we offer possible remedies. As for the slow convergence, symmetric Nitsche seems comparable or even slightly superior in performance. Furthermore, there are multiple implementations of findinf a normal vectors. See the .pdf and references listed there for details. In short, one can measure the normal at the element sides or vertices. Both choices have some advantages, but especially in the case of symmetric Nitsche, the "vertex normal" might be preferable. Our implementation uses the "side normal", so symmetric Nitsche with "vertex normal" might speed-up the computations. More suggestions and comments are listed in the .pdf document.

As for the seeming convergence. If this happens -- depends on material parameters, phonon flow momentum, stabilisation parameter, mesh size and (very importantly) geometry -- it seems to enforce near no-slip boundary conditions. Therefore, if one refines the mesh and observes no change in the solutions and a nontrivial flow along the boundary, the solution likely is the correct one. We, however, highly recommend the check of slipping on the boundary. Finally, it is noteworthy that only a change in geometry can have a massive impact. In our computations, Nitsche method gives perfect results and fast convergence rate if applied to the straight walls of the channel (which are tangential to the flow), while it gives the worryingly bad results if applied to the cylindrical obstacle, which is much curved and much more exposed to the flow.

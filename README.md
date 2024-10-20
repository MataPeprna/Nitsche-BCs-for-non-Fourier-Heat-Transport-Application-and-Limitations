# Nitsche-Boundary-Conditions-for-non-Fourier-Heat-Transport
Non-Fourier heat conduction is a field of complex models describing heat transport in very specific a challenging settings (i.e. extremely small sizes etc.). This project presents a way to model slip boundary conditions for a complex non-Fourier heat conduction model and implements it numerically.

__________________

Fourier law states (roughly speaking) that heat travels in the direction, where temperature is decreasing the fastest. This sounds logical and it is usually true. However, in certain complex settings like very small size scale etc., this law is not valid and heat conduction rather resembles flow of a fluid (for instance wind). This sort of heat transport can be described by a rather complex field called phonon hydrodynamics. (Phonons are quasi particles describing for example heat in a crystal and hydrodynamics is a science describing flows od fluids.)

One of the recently developed models of phonon hydrodynamics is particularly similar to the famous Navier-Stokes equations. This projects implements full-slip* (and also Navier-slip) boundary conditions to the model by the so-called Nitsche method and studies, how to choose a stabilisation parameter to optimise performance.

This project contains a .pdf summary of the results which serve as recommendations on how to choose the parameter and working codes that give the results from the .pdf. In order to make the codes work, one needs to download legacy FEniCS 2019.1.0, see [this site](https://fenicsproject.org/download/archive/) for download.

*Full-slip says that there is no flow through the boundary, but the boundary is perfectly slippery and does not slow down the flowin the tangential direction (along itself).

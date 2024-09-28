# ThermalModelElOro.jl

This folder contains the code used to produce the 2 thermal models of the El Oro case study, a convection driven model and an advection model. The code is organised as follows: the folder src contains the code of the unregistered package ThermalModelElOro. This package uses the package DifferentialEquations.jl to model a 1D thermal profile. 

The file 1D_diffusion_advection_convection.jl contains the code to run both models. It will generate a figure similar to Figure 7.

To run the model, you can use the following command, if you are in the ThermalModelElOro.jl folder:

```
julia --project 1D_diffusion_advection_convection.jl
```
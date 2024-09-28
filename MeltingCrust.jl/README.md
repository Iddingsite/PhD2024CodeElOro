# MeltingCrust.jl

This folder contains the code used to generate the results of the coupling between MAGEMin and the convection thermal model. The code is organised as follows: the folder src contains the code of the unregistered package MeltingCrust. This package uses the package Differentialequations.jl at its core to model a 1D thermal profile. Callback functions have been defined to be able to interact with the thermal model at specific timesteps (see DiffEqCallbacks.jl package for more details). This includes calling MAGEMin to calculate the phase equilibrium of each node of the model.

The folder ElOro contains the initial conditions and the files used to run the models for the El Oro case study. 

This can be run using the following command, if you are in the MeltingCrust.jl folder:

```
julia --project ElOro/El_Oro_convection_model.jl
```

This will run a model with the parameters used for the Biotite model. All the results will be saved in a hdf5 file in the folder ElOro. This can be read with HDFview or any programming language. This contains all the timesteps of the model, with the composition of melt extracted, the temperature, phase stables in the 1D model, etc..

To change the parameters of the model, the line number 22 of the file El_Oro_convection_model.jl can be modified:

```julia
# this is the important line where we define the parameters for the MAGEMin interaction
magemin_interact = MAGEMinInteraction(grid=grid_model, crystal_fractionation=false, frac_crystallisation=[0.0, 0.0, 0.0], flux_melting=false, h2o_keep= 0.00, ms_melt_only=false, grt_frac_remove=0.0, crd_frac_remove = 0.0)
```
For instance, to run a model with only extraction when muscovite is stable and to produce fractionation of 10% of plagioclase, the line should be modified like this:

```julia
# this is the important line where we define the parameters for the MAGEMin interaction
magemin_interact = MAGEMinInteraction(grid=grid_model, crystal_fractionation=true, frac_crystallisation=[0.1, 0.0, 0.0], flux_melting=false, h2o_keep= 0.00, ms_melt_only=true, grt_frac_remove=0.0, crd_frac_remove = 0.0)
```
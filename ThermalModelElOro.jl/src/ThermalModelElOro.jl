module ThermalModelElOro


using Reexport
@reexport using OrdinaryDiffEq, DiffEqCallbacks, LinearSolve
using ProgressBars
using Plots
using GeoParams
@reexport using Parameters
@reexport using Unitful
using Logging: global_logger
using TerminalLoggers: TerminalLogger



function __init__()
    # initialise global logger for DifferentialEquations
    global_logger(TerminalLogger())
end

include("input/InitialConditions.jl")
include("thermalproperties.jl")
include("semidiscretization.jl")
include("callbacks/boundarychange.jl")
include("callbacks/callthermo.jl")

export CreateGrid, ThermalVariables, PhysicalProperties
export thermal_diffusivity!, heat_capacity!, thermal_conductivity!, convective_heat_flux!, convective_thermal_conductivity!
export semi_discretisation
export change_boundary!, change_boundary_convection_flux!, change_boundary_convection_gabbro!, change_boundary_convection_blueschist!
export remove_and_emplace_melt_func

end # module ThermalModelElOro

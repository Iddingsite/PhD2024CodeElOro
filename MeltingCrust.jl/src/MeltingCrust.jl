module MeltingCrust

using Reexport

@reexport using MAGEMin_C
@reexport using OrdinaryDiffEq, DiffEqCallbacks, LinearSolve
using ProgressBars
using Plots
@reexport using Parameters
@reexport using Unitful
using Logging: global_logger
using TerminalLoggers: TerminalLogger
using GeoParams


# global data_MAGEMin = Initialize_MAGEMin("mp", verbose=false);
# export data_MAGEMin

function __init__()
    # initialise global logger for DifferentialEquations
    global_logger(TerminalLogger())
end

include("utils/convert_to_wt.jl")
include("input/InitialConditions.jl")
include("thermalproperties.jl")
include("semidiscretization.jl")
include("MAGEMincall/call_MAGEMin.jl")
include("callbacks/updateproperties.jl")
include("callbacks/boundarychange.jl")
include("callbacks/callthermo.jl")
include("callbacks/output.jl")

export CreateGrid, ThermalVariables, PhysicalProperties, MAGEMinInteraction, Compositions
export thermal_diffusivity!, heat_capacity!, thermal_conductivity!, convective_heat_flux!, convective_thermal_conductivity!
export semi_dicretisation
export convert_compo_to_wt_pct!, convert_compo_to_wt_dry_pct!
export change_boundary!, thermal_param_call_func, minimization_calls_magemin!, change_boundary_convection_gabbro!, change_boundary_convection_gabbro!, change_boundary_convection_blueschist!, change_boundary_convection_flux!
export Initialize_MAGEMin, find_ox_index!, extract_H2O_MAGEMin!, extract_melt_MAGEMin!, update_compo_MAGEMin_i!
export stencil_diffusion!
export save_data

end # module MeltingCrust

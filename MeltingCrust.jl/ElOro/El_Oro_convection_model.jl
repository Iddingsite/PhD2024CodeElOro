using MeltingCrust
using DelimitedFiles
using Plots
using GeoParams
using Dates


println("Thread number: $(Threads.nthreads())")

cd(@__DIR__)

include("initial_conditions_El_Oro.jl")

# physical properties
grid_model = CreateGrid(nz=200, Lz=28u"km", tfinal=10.5u"Myr")

thermal_param = ThermalVariables(grid=grid_model,heat_flux_moho=0, time_cooling=10 * 1e6 * 365.25 * 24 * 3600, convection=[false])
physical_prop = PhysicalProperties(grid=grid_model)
compositions = Compositions(grid=grid_model)

# this is the important line where we define the parameters for the MAGEMin interaction
magemin_interact = MAGEMinInteraction(grid=grid_model, crystal_fractionation=false, frac_crystallisation=[0.0, 0.0, 0.0], flux_melting=false, h2o_keep= 0.00, ms_melt_only=false, grt_frac_remove=0.0, crd_frac_remove = 0.0)

# print the parameters for the MAGEMin interaction
println("Retains the melt in the system: $(magemin_interact.ms_melt_only)")
println("Fraction of crystallisation: $(magemin_interact.frac_crystallisation)")
println("Fraction of garnet removed: $(magemin_interact.grt_frac_remove)")
println("Fraction of cordierite removed: $(magemin_interact.crd_frac_remove)")
println("Water retained in the system: $(magemin_interact.h2o_keep)")


@unpack Δz, Δz_, zc, z, nz, depthc = grid_model
@unpack T0, α, Cp, k, q_conv, k_conv, ρ, Ra, Nu = thermal_param
@unpack Plith = physical_prop


initial_T!(T0, zc);
initial_Plith!(Plith, ρ, z)

# add 1 kbar of pressure to prevent occurence of Crd
Plith .+= 1e8

# define material parameters using GeoParams
MatParam = (SetMaterialParams(;Name="Crust", Phase=1,
                        Density   = ConstantDensity(ρ=2700kg/m^3),
                        HeatCapacity = T_HeatCapacity_Whittington(),
                        Conductivity = T_Conductivity_Whittington_parameterised()),)

# define phase array with 1 for z bellow 25_000 and 2 above
phases = ones(size(T0,1),1)
phases .= 1

# define density
compute_density!(ρ, MatParam, phases, (; T=T0))
initial_Nu!(Nu, Ra, zc)

compute_heatcapacity!(Cp, MatParam, phases, (;T=T0))
compute_conductivity!(k, MatParam, phases, (;T=T0))

convective_heat_flux!(q_conv, Nu, k)
convective_thermal_conductivity!(k_conv, q_conv, k)

# initial compositions
@unpack compo_solid, compo_solid_matrix, compo_pelite_forshaw, compo_pelite_forshaw_ox, compo_arkose_pettijohn,  compo_arkose_pettijohn_ox, compo_quartzite_clarke, compo_quartzite_clarke_ox, compo_pelite_riel,  compo_pelite_riel_ox, compo_semipelite_riel, compo_semipelite_riel_ox, compo_psammite_riel, compo_psammite_riel_ox, compo_solid_name, compo_pelite_dominguez_ox, compo_semipelite_dominguez_EOET_8_ox = compositions

#! mole of oxides for MAGEMin

layers_compo!(compo_solid, compo_solid_name, compo_pelite_dominguez_ox, compo_semipelite_dominguez_EOET_8_ox, compo_quartzite_clarke_ox)

@unpack depth = grid_model
@unpack depth_min_thermo, depth_max_thermo, miniz_calls, ΔT,oxides_name = magemin_interact
@unpack Plith = physical_prop

data_MAGEMin = Initialize_MAGEMin("mp", verbose=false, solver=1);

println("Initialising the initial conditions with MAGEMin")
initial_call_MAGEMin!(compo_solid, grid_model, physical_prop, thermal_param, magemin_interact, T0, data_MAGEMin)

datetoday = string(Dates.today())
path_save = "ElOro_frac_$(magemin_interact.frac_crystallisation)_flux_melting_$(magemin_interact.flux_melting)_water_keep_$(magemin_interact.h2o_keep)_ms_melt_$(magemin_interact.ms_melt_only)_compo_dominguez_grt_remove_$(magemin_interact.grt_frac_remove)_crd_remove_$(magemin_interact.crd_frac_remove)_$(datetoday).h5"
# path_save = "melting_dominguez/test.h5"

p = (grid=grid_model, thermal_param=thermal_param, compositions=compositions, physical_prop=physical_prop, magemin_interact=magemin_interact, data_MAGEMin=data_MAGEMin, path_save=path_save, geoparam=(MatParam, phases))
t = (0.0, grid_model.tfinal)


ΔT_step = ustrip.(u"s", vcat(
    (0u"Myr":1000u"yr":8u"Myr"),
    (8u"Myr"+100u"yr":100u"yr":8.5u"Myr"),
    (8.5u"Myr"+250u"yr":250u"yr":10u"Myr"),
    (10u"Myr"+1000u"yr":1000u"yr":12u"Myr"),
    (12u"Myr"+500000u"yr":500000u"yr":20u"Myr"),
))

save_step = ustrip.(u"s", vcat(
    (0u"Myr":100000u"yr":8u"Myr"),
    (8u"Myr"+500u"yr":500u"yr":12u"Myr"),
    (12u"Myr"+100000u"yr":100000u"yr":20u"Myr")
))

emplacement_gabbro = PresetTimeCallback([8 * 1e6 * 365.25 * 24 * 3600],change_boundary_convection_gabbro!, save_positions=(false,false))
cooling_gabbro = PresetTimeCallback([thermal_param.time_cooling],change_boundary_convection_blueschist!, save_positions=(false,false))
flux_boundary_call = FunctionCallingCallback(change_boundary_convection_flux!; funcat=Vector{Float64}(), func_everystep=true, func_start = false, tdir=1);

miniz_nb_check = PresetTimeCallback(ΔT_step,minimization_calls_magemin!, save_positions=(false,false))
save_data_callback = PresetTimeCallback(save_step, save_data, save_positions=(false,false))

callbacks = CallbackSet(emplacement_gabbro, cooling_gabbro, flux_boundary_call, miniz_nb_check, save_data_callback)

# ROCK2 using DiffEq
prob = ODEProblem(semi_dicretisation, T0, t, p);
@time sol = solve(prob, Tsit5(), progress = true, progress_steps = 1, save_start=true, callback= callbacks, save_everystep=true, abstol=1e-6,reltol=1e-6);

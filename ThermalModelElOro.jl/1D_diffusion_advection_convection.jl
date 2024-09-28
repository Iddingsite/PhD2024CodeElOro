using ThermalModelElOro
using GeoParams
using OrdinaryDiffEq
using Plots
using BenchmarkTools
using ProgressBars

grid_model = CreateGrid(nz=200, Lz=30u"km", tfinal=5u"Myr")

thermal_param = ThermalVariables(grid=grid_model, M=10,heat_flux_moho=-150*1e-3)
physical_param = PhysicalProperties(grid=grid_model)

@unpack Δz, Δz_, zc, depthc, nz, tfinal = grid_model
@unpack T0, α, Cp, k, q_conv, k_conv, ρ, Ql, A, Nu, Ra, input_gabbro, crust_array, convection, advection = thermal_param
@unpack Plith, ϕ, dϕdT = physical_param

convection .= false
advection .= false
input_gabbro .= false

# Initial Temperature in K
for i in axes(zc,1)
    T0[end - i + 1] = 0.025 * zc[i] + 0 + 273.15
end

T0_adv = copy(T0)

# rhyolite melt curve
p_rhyolite = MeltingParam_Quadratic(T_s=(680 + 273.15), T_l=(1150 + 273.15))
p_basalt = MeltingParam_Quadratic(T_s=(820 + 273.15), T_l=(1200 + 273.15))


# define material parameters using GeoParams
MatParam = (SetMaterialParams(;Name="FelsicCrust", Phase=1,
                        Density   = ConstantDensity(ρ=2700kg/m^3),
                        HeatCapacity = T_HeatCapacity_Whittington(),
                        Conductivity = T_Conductivity_Whittington_parameterised()),
            SetMaterialParams(;Name="MaficCrust", Phase=2,
                        Density   = ConstantDensity(ρ=3300kg/m^3),
                        HeatCapacity = T_HeatCapacity_Whittington(),
                        Conductivity = T_Conductivity_Whittington_parameterised()))

# define phase array with 1 for z bellow 30_000 and 2 above
phases = ones(size(T0,1),2)
for i in axes(depthc,1)
    if depthc[i] >= 25_000
        phases[i,1] = 0
        phases[i,2] = 1
    else
        phases[i,2] = 0
    end
end

# define density
compute_density!(ρ, MatParam, phases, (; T=T0))

function initial_Plith!(Plith, ρ, Δz)
    # Initial Plith
    for i in axes(Plith,1)
        if i == 1
            Plith[i] = 0
        else
            Plith[i] = Plith[i-1] + ρ[i] * 9.81 * Δz
        end
    end
end

initial_Plith!(Plith, reverse(ρ), Δz)

compute_heatcapacity!(Cp, MatParam, phases, (; T=T0))
# define thermal conductivity

compute_conductivity!(k, MatParam, phases, (; T=T0))


# define a named tuple containing the parameters for the ODE solver
p = (geoparam=(MatParam, phases), grid=grid_model, thermal_param=thermal_param, physical_param=physical_param)
# define the initial and final time
t = (0, tfinal)

flux_boundary_call = PresetTimeCallback([thermal_param.time_cooling],change_boundary!)

# define the ODE problem
prob = ODEProblem(semi_discretisation, T0, t, p);

melt_emplacement = FunctionCallingCallback(remove_and_emplace_melt_func; funcat=Vector{Float64}(), func_everystep=true, func_start = false, tdir=1);

callbacks = CallbackSet(melt_emplacement)

# solve with ROCK2 (stabilized explicit solver) using DiffEq
@time sol = solve(prob, ROCK2(), progress = true, progress_steps = 1, save_start=true, abstol=1e-6,reltol=1e-6, callback=callbacks);

melt_emplacement = copy(thermal_param.melt_emplaced)


# from Riel et al. 2016
PT_points_T = [799.4632313472895,809.9469940955452,815.1888754696728,878.0914519592058,880.7123926462697,757.5281803542675,707.7303073000539,600.2717391304349,583.2356246645198]
PT_points_depth = [17.242374103771585,18.87888586610666,20.26681356327691,21.571879905392223,23.47769107165585,16.84878266726062,15.502285647617837,12.581422881931188,11.607801960035639]
P_uncertainties = [0.6,0.6,0.6,0.6,0.6,0.6,2.5,2.5,2.5,]
T_uncertainties = [55,55,55,55,55,55,55,55,55]


plot(grid_model.depthc.*1e-3, sol[end] .- 273.15, label = "Initial conditions", dpi=400, title = "Initial conditions")
scatter!(PT_points_depth, PT_points_T, yerror=T_uncertainties, xerror=P_uncertainties, label = "PT points", marker = :circle, markersize = 3, legend = :topleft)

sol_next = deepcopy(sol[end])

# find index of 25_000
boundary_top = findall(grid_model.depthc .>= 25_000)[end]
boundary = findall(grid_model.depthc .>= 28_000)[end]

sol_next[boundary:boundary_top] .= 1250 + 273.15

sol_next = copy(sol_next[boundary:end])

grid_model_2 = CreateGrid(nz=size(sol_next,1), Lz=28u"km", tfinal=10u"Myr")
thermal_param = ThermalVariables(grid=grid_model_2, M=1, time_cooling=0.3*(3600*24*365.25*1e6), heat_flux_moho=-150*1e-3, neumann_bottom=false)
physical_param = PhysicalProperties(grid=grid_model_2)

@unpack Δz, Δz_, zc, depthc, nz, tfinal = grid_model_2
@unpack α, Cp, k, q_conv, k_conv, ρ, Ql, A, Nu, Ra, input_gabbro, crust_array, convection, advection = thermal_param
@unpack Plith, ϕ, dϕdT = physical_param


input_gabbro .= false
convection .= false
advection .= false

# define phase array with 1 for z bellow 30_000 and 2 above
phases = ones(size(sol_next,1),2)
for i in axes(depthc,1)
    if depthc[i] >= 25_000
        phases[i,1] = 0
        phases[i,2] = 1
    else
        phases[i,2] = 0
    end
end

# define density
compute_density!(ρ, MatParam, phases, (; T=sol_next))


# define a named tuple containing the parameters for the ODE solver
p = (geoparam=(MatParam, phases), grid=grid_model_2, thermal_param=thermal_param, physical_param=physical_param)
# define the initial and final time
t = (0, tfinal)

flux_boundary_call = PresetTimeCallback([thermal_param.time_cooling],change_boundary!)

# define the ODE problem
prob = ODEProblem(semi_discretisation, sol_next, t, p);

@time sol_2 = solve(prob, ROCK2(), progress = true, progress_steps = 1, save_start=true, abstol=1e-6,reltol=1e-6, callback=flux_boundary_call);



plot(grid_model.depthc.*1e-3, sol_2[end] .- 273.15, label = "Initial conditions", dpi=400, title = "Initial conditions")
scatter!(PT_points_depth, PT_points_T, yerror=T_uncertainties, xerror=P_uncertainties, label = "PT points", marker = :circle, markersize = 3, legend = :topleft)


grid_model_conv = CreateGrid(nz=200, Lz=28u"km", tfinal=20u"Myr")

thermal_param_conv = ThermalVariables(grid=grid_model_conv, heat_flux_moho=-150*1e-3, time_cooling=10 * 1e6 * 365.25 * 24 * 3600)
physical_param_conv = PhysicalProperties(grid=grid_model_conv)


@unpack Δz, Δz_, zc, nz, tfinal, depthc = grid_model_conv
@unpack T0, α, Cp, k, q_conv, k_conv, ρ, A, Nu, Ra, input_gabbro, crust_array, convection, advection = thermal_param_conv
@unpack Plith, ϕ = physical_param_conv

convection .= false
advection .= false
input_gabbro .= false

# Initial Temperature in K
for i in axes(zc,1)
    T0[i] = 0.025*zc[i] + 0 + 273.15
end
# T0[end] = 1250 .+ 273.15
T0 .= reverse(T0, dims=1)

# define material parameters using GeoParams
MatParam = (SetMaterialParams(;Name="Crust", Phase=1,
                        Density   = ConstantDensity(ρ=2700kg/m^3),
                        HeatCapacity = T_HeatCapacity_Whittington(),
                        Conductivity = T_Conductivity_Whittington_parameterised()),
            SetMaterialParams(;Name="Mantle", Phase=2,
                        Density   = ConstantDensity(ρ=3300kg/m^3),
                        HeatCapacity = T_HeatCapacity_Whittington(),
                        Conductivity = T_Conductivity_Whittington_parameterised()))

# define phase array with 1 for z bellow 25_000 and 2 above
phases = ones(size(T0,1),2)
for i in axes(zc,1)
    if reverse(zc)[i] >= 25_000
        phases[i,1] = 0
        phases[i,2] = 1
    else
        phases[i,2] = 0
    end
end

# define density
compute_density!(ρ, MatParam, phases, (; T=T0))

function initial_Plith!(Plith, ρ, Δz)
    # Initial Plith
    for i in axes(Plith,1)
        if i == 1
            Plith[i] = 0
        else
            Plith[i] = Plith[i-1] + ρ[i] * 9.81 * Δz
        end
    end
end

initial_Plith!(Plith, reverse(ρ), Δz)

compute_heatcapacity!(Cp, MatParam, phases, (; T=T0))
# define thermal conductivity

compute_conductivity!(k, MatParam, phases, (; T=T0))

for i in axes(zc,1)
    if reverse(zc)[i] >= 17_000 && reverse(zc)[i] <= 25_000
        Nu[i] = 0.120*Ra^(1/3)
    end
end

convective_heat_flux!(q_conv, Nu, k)
convective_thermal_conductivity!(k_conv, q_conv, k)

# define a named tuple containing the parameters for the ODE solver
p = (geoparam=(MatParam, phases), grid=grid_model_conv, thermal_param=thermal_param_conv, physical_param=physical_param_conv)
# define the initial and final time
t = (0, tfinal)

# define callback
emplacement_gabbro = PresetTimeCallback([8 * 1e6 * 365.25 * 24 * 3600],change_boundary_convection_gabbro!)

cooling_gabbro = PresetTimeCallback([thermal_param_conv.time_cooling],change_boundary_convection_blueschist!)


# define the ODE problem
prob = ODEProblem(semi_discretisation, T0, t, p);

# solve with ROCK2 (stabilized explicit solver) using DiffEq

# flux_boundary_call = PresetTimeCallback([thermal_param.time_cooling],change_boundary!)
flux_boundary_call = FunctionCallingCallback(change_boundary_convection_flux!; funcat=Vector{Float64}(), func_everystep=true, func_start = false, tdir=1);
callbacks = CallbackSet(emplacement_gabbro, cooling_gabbro, flux_boundary_call)

@time sol_conv = solve(prob, ROCK2(), progress = true, progress_steps = 1, save_start=true, abstol=1e-6,reltol=1e-6, callback= callbacks);



plot(grid_model_2.depthc.*1e-3, sol_2[end] .- 273.15, label = "Initial conditions", dpi=400, title = "Initial conditions")
plot!(grid_model.depthc .* 1e-3, sol[end] .- 273.15, label = "Advection")


# from Riel et al. 2016
PT_points_T = [799.4632313472895,809.9469940955452,815.1888754696728,878.0914519592058,880.7123926462697,757.5281803542675,707.7303073000539,600.2717391304349,583.2356246645198]
PT_points_depth = [17.242374103771585,18.87888586610666,20.26681356327691,21.571879905392223,23.47769107165585,16.84878266726062,15.502285647617837,12.581422881931188,11.607801960035639]
P_uncertainties = [0.6,0.6,0.6,0.6,0.6,0.6,2.5,2.5,2.5,]
T_uncertainties = [55,55,55,55,55,55,55,55,55]

# create a figure with 3 subplots
layout = @layout [a b; c]

maximum_T_conv_1 = zeros(size(T0,1))
maximum_T_conv_2 = zeros(size(sol_next,1))

# loop over sol and sol_next to get the maximum values and feed it to maximum_T_conv_1 and maximum_T_conv_2

for t = sol.t
    for i = axes(maximum_T_conv_1,1)
        maximum_T_conv_1[i] = max(maximum_T_conv_1[i],sol(t)[i])
    end
end

maximum_T_conv_2 .= maximum_T_conv_1[boundary:end]

for t = sol_2.t
    for i = axes(maximum_T_conv_2,1)
        maximum_T_conv_2[i] = max(maximum_T_conv_2[i],sol_2(t)[i])
    end
end



# plot the temperature evolution of advection in a, convection in b and the points at different depths in c

using ColorSchemes


# create a colorbar with 30 values from 0 to grid_model.tfinal+grid_model_2.tfinal
colors_prograde = [get((ColorSchemes.OrRd_4), i/11) for i in 0:10]

colors_retrograde = [get((ColorSchemes.Blues_3), i/7) for i in 0:6]

# plot condition initial
a = plot(sol[1] .- 273.15, grid_model.depthc.*1e-3,  label = "Initial conditions", dpi=400, title = "Total duration: $(Int(round(((grid_model.tfinal+grid_model_2.tfinal) / (1e6*3600*365.25*24));digits=0))) Myr\nTotal intrusion thickness: $(round(melt_emplacement*1e-3;digits=2)) km, M=10", yflip=true, xlabel="Temperature [°C]", ylabel="Depth [km]", color=:black, legend=:topright, linewidth=1.2, titlefontsize=6, ylims=(0,30), legendfont=font(6), xlims=(0,1260), minorgrid = false, xminorticks=2, yminorticks=5, background_color_inside=:gray98, gridstyle=:dash) 

step_s = 0.5*3600*24*365.25*1e6
tfinal_plot = grid_model.tfinal+thermal_param.time_cooling
tfinal_plot_Ma = tfinal_plot/(3600*24*365.25*1e6)

range_plot = range(start=0, stop=tfinal_plot, step=step_s )

for (i,t) = tqdm(enumerate(range_plot))

    if t < grid_model.tfinal
        plot!(a, sol(t)[:] .- 273.15, reverse(grid_model.zc)*.001, label = "", alpha=1, linewidth=1.2, color=colors_prograde[i])
    else
        plot!(a, sol_2(t-grid_model.tfinal)[:] .- 273.15, reverse(grid_model_2.zc)*.001, label = "", alpha=1, linewidth=1.2, color=colors_prograde[i])
    end
end


tfinal_plot = grid_model.tfinal+grid_model_2.tfinal
tfinal_plot_Ma = tfinal_plot/(3600*24*365.25*1e6)
step_s = 1.5*3600*24*365.25*1e6

range_plot = range(start=grid_model.tfinal+thermal_param.time_cooling, stop=tfinal_plot, step=step_s )

for (i,t) = tqdm(enumerate(range_plot))
    plot!(a, sol_2(t-grid_model.tfinal)[:] .- 273.15, reverse(grid_model_2.zc)*.001, label = "", alpha=1, linewidth=1.4, color=colors_retrograde[i], linestyle=:dot)
end

# plot maximum maximum_T_conv_2
scatter!(a, PT_points_T, PT_points_depth, label = "P–T estimates from\nRiel et al., 2013, 2016", marker = :pentagon, markersize = 3, yflip = true, yerror=P_uncertainties, xerror=T_uncertainties, color="gold2", markerstrokewidth=.5)
plot!(a, maximum_T_conv_2[:] .- 273.15, grid_model_2.depthc.*1e-3, label = "Maximum T", color=:firebrick, linewidth=1.8)
# plot condition initial
plot!(a,sol[1] .- 273.15, grid_model.depthc.*1e-3,  label = "",color=:black, legend=:topright, linewidth=1.8)

# plot initial conditions

current()

# do the same thing for convection

# create a colorbar with 30 values from 0 to grid_model.tfinal+grid_model_2.tfinal
colors_prograde = [get((ColorSchemes.OrRd_4), i/22) for i in 0:21]

colors_retrograde = [get((ColorSchemes.Blues_3), i/6) for i in 0:5]


b = plot(sol_conv[1] .- 273.15, grid_model_conv.depthc.*1e-3,  label = "Initial conditions", dpi=400, title = "Total duration: $(Int(round(grid_model_conv.tfinal / (1e6*3600*365.25*24);digits=0))) Myr\nDuration convection: 2 Myr", yflip=true, xlabel="Temperature [°C]", ylabel="", color=:black, legend=:topright, linewidth=1.2, titlefontsize=6, ylims=(0,30), legendfont=font(6),ytickfontcolor = RGBA(0,0,0,0), xlims=(0,1260), minorgrid = false, xminorticks=2, yminorticks=5, background_color_inside=:gray98, gridstyle=:dash)

step_s = 0.5*3600*24*365.25*1e6
tfinal_plot = grid_model_conv.tfinal
tfinal_plot_Ma = tfinal_plot/(3600*24*365.25*1e6)

range_plot = range(start=0, stop=thermal_param_conv.time_cooling, step=step_s )

for (i,t) = tqdm(enumerate(range_plot))

    plot!(b, sol_conv(t)[:] .- 273.15, reverse(grid_model_conv.zc)*.001, label = "", alpha=1, linewidth=1.2, color=colors_prograde[i])
end

step_s = 2*3600*24*365.25*1e6
range_plot = range(start=thermal_param_conv.time_cooling, stop=tfinal_plot, step=step_s )

for (i,t) = tqdm(enumerate(range_plot))
    plot!(b, sol_conv(t)[:] .- 273.15, reverse(grid_model_conv.zc)*.001, label = "", alpha=1, linewidth=1.4, color=colors_retrograde[i], linestyle=:dot)
end

maximum_T = zeros(size(T0,1))

for t = sol_conv.t
    for i = axes(maximum_T,1)
        maximum_T[i] = max(maximum_T[i],sol_conv(t)[i])
    end
end

scatter!(b, PT_points_T, PT_points_depth, label = "P–T estimates from\nRiel et al., 2013, 2016", marker = :pentagon, markersize = 3, yflip = true, yerror=P_uncertainties, xerror=T_uncertainties, color="gold2", markerstrokewidth=.5)
plot!(b, maximum_T[:] .- 273.15, grid_model_conv.depthc.*1e-3, label = "Maximum T", color=:firebrick, linewidth=1.8)
plot!(b, sol_conv[1] .- 273.15, grid_model_conv.depthc.*1e-3,  label = "",color=:black, legend=:topright, linewidth=1.8)


current()

depth_plot = [14000, 17000, 19000, 22000]

# plot each depth in a tvsT plot
t_total_adv = vcat(sol.t, sol_2.t .+ grid_model.tfinal) ./ (3600*24*365.25*1e6)
t_total_conv = sol_conv.t ./ (3600*24*365.25*1e6)

# find maximum temperature and time at each depth 
point_T_adv_max = [zeros(size(depth_plot,1)), zeros(size(depth_plot,1))]
point_T_conv_max = [zeros(size(depth_plot,1)), zeros(size(depth_plot,1))]

values_T_adv = zeros(size(sol.u,1) + size(sol_2.u,1), size(depth_plot,1))
values_T_conv = zeros(size(sol_conv.u,1), size(depth_plot,1))

# find the closest depth to the depth_plot
index_sol = zeros(Int, size(depth_plot,1))
index_sol_2 = zeros(Int, size(depth_plot,1))
index_sol_conv = zeros(Int, size(depth_plot,1))

for (j,depth) = enumerate(depth_plot)
    index_emplacement = argmin(abs.(grid_model.depthc .- (depth)))
    index_sol[j] = index_emplacement
end

# now for index_sol_2
for (j,depth) = enumerate(depth_plot)
    index_emplacement = argmin(abs.(grid_model_2.depthc .- (depth)))
    index_sol_2[j] = index_emplacement
end

for (j,depth) = enumerate(depth_plot)
    index_emplacement = argmin(abs.(grid_model_conv.depthc .- (depth)))
    index_sol_conv[j] = index_emplacement
end


# iterate through values_T_adv
for i in axes(values_T_adv,2)
    values_T_adv[:,i] = [getindex.(sol.u,index_sol[i]); getindex.(sol_2.u,index_sol_2[i])]
end

for i in axes(values_T_conv,2)
    values_T_conv[:,i] = getindex.(sol_conv.u,index_sol_conv[i])
end

# find maximum temperature and time at each depth in point_T_adv_max and point_T_conv_max using argmax
for (i,depth) = enumerate(depth_plot)
    point_T_adv_max[1][i] = maximum(values_T_adv[:,i])
    point_T_adv_max[2][i] = t_total_adv[argmax(values_T_adv[:,i])]
end

for (i,depth) = enumerate(depth_plot)
    point_T_conv_max[1][i] = maximum(values_T_conv[:,i])
    point_T_conv_max[2][i] = t_total_conv[argmax(values_T_conv[:,i])]
end




c = plot(title="", xlabel="Time [Myr]", ylabel="Temperature [°C]", dpi=400, size=(800,500), left_margin = [3Plots.mm 0Plots.mm], minorgrid = false, xminorticks=5, yminorticks=2, background_color_inside=:gray98, xlims=(0,20), ylims=(200,900), gridstyle=:dash)

color_c_adv = [get((ColorSchemes.matter), i/(4)) for i in 0:3]
color_c_conv = reverse([get(((ColorSchemes.davos10)), i/(4)) for i in 0:3])

for (i,depth) = enumerate(depth_plot)
    plot!(c, t_total_adv .+ 4.7, values_T_adv[:,i] .- 273.15, label = "$(depth*1e-3) km", color=color_c_adv[i], linewidth=2)
end

for (i,depth) = enumerate(depth_plot)
    plot!(c, t_total_conv, values_T_conv[:,i] .- 273.15, label = "$(depth*1e-3) km", color=color_c_conv[i], linewidth=2, linestyle=:dashdot)
end

# for values_T_adv, replace the values at 17000 and 19000 by the temperature at 10 Ma from sol
point_T_adv_max[1][2] = sol_2[1][findall(grid_model_2.depthc .>= 17000)[end]]
point_T_adv_max[2][2] = ustrip(5u"Myr")
point_T_adv_max[1][3] = sol[end][findall(grid_model.depthc .>= 19000)[end]]
point_T_adv_max[2][3] = ustrip(5u"Myr")

# plot maximum T and time for each depth with a star
scatter!(c, point_T_adv_max[2] .+ 4.7, point_T_adv_max[1] .- 273.15, label = "", marker = :star5, markersize = 6, color=color_c_adv)
scatter!(c, point_T_conv_max[2], point_T_conv_max[1] .- 273.15, label = "", marker = :star5, markersize = 6, color=color_c_conv)


current()

# plot all the plots in a single figure
plot(a,b,c, layout=layout, dpi=400, size=(600,800), left_margin = [3Plots.mm 0Plots.mm])

savefig("El_Oro_Thermal_models.pdf")

using HDF5

const yr_to_myr = 365.25 * 24 * 3600 * 1e6

"""
    set_attributes_and_groups(t, tcurrent, Δt, magemin_interact, integrator, Plith, compo_solid_matrix)

This function sets attributes and groups for the given HDF5 group `t`. It takes the current time `tcurrent`, time step `Δt`, `magemin_interact` object, `integrator` object, lithostatic pressure `Plith`, and `compo_solid_matrix` as arguments. It creates groups and datasets in the HDF5 file and sets their attributes.
"""
function set_attributes_and_groups(t, tcurrent, Δt, magemin_interact, integrator, Plith, compo_solid_matrix, compo_solid_name)
    @unpack compo_rock_wt, compo_i_g, total_volume_melt, compo_melt_wt_all, compo_melt_wt_all_index, compo_melt_all_time, compo_melt_wt_dry_all, compo_assemblage, compo_melt_all_volume, compo_melt_all_unit = magemin_interact

    # check if compo_melt_wt_all is not empty
    if !isempty(compo_melt_wt_all)
        compo_melt_wt_save = hcat(compo_melt_wt_all...)'
        compo_melt_wt_save_dry = hcat(compo_melt_wt_dry_all...)'
    end


    tcurrent /= yr_to_myr
    Δt /= yr_to_myr
    compo_melt_all_time /= yr_to_myr

    attributes(t)["Time(Myr)"] = tcurrent
    attributes(t)["CurrentDt(Myr)"] = Δt
    attributes(t)["TotalMinimisationCall"] = magemin_interact.miniz_calls
    attributes(t)["MeltTotalVolume(m³)"] = total_volume_melt


    # Create groups and set attributes for each group
    groups = ("T(°C)", "Plith(kbar)", "Rock", "Melt")
    for group in groups
        grp = create_group(t, group)
        attributes(grp)["DataType"] = "Scalar"
        attributes(grp)["Center"] = "Node"
    end

    # Create datasets and set attributes for each dataset
    t["T(°C)"]["T(°C)"] = convert(Array{Float32}, integrator.u .- 273.15)
    t["Plith(kbar)"]["Plith(kbar)"] = convert(Array{Float32}, Plith ./ 1e8)
    t["Rock"]["CompositionRock(mol%)"] = convert(Array{Float32}, compo_solid_matrix)
    t["Rock"]["CompositionRock(wt%)"] = convert(Array{Float32}, compo_rock_wt)
    t["Rock"]["CompositionAssemblage(mol%ofOxides)"] = convert(Array{Float32}, compo_assemblage')
    t["Rock"]["UnitName"] = compo_solid_name
    t["Rock"]["ShrinkingVolume(frac)"] = magemin_interact.vol_frac_system_shrinking_tot
    if !isempty(compo_melt_wt_all)
        t["Melt"]["CompositionMelt(wt%)"] = convert(Array{Float32}, compo_melt_wt_save)
        t["Melt"]["CompositionMeltDry(wt%)"] = convert(Array{Float32}, compo_melt_wt_save_dry)
        t["Melt"]["IndexGridMelt"] = compo_melt_wt_all_index
        t["Melt"]["TimeMelt(Myr)"] = compo_melt_all_time
        t["Melt"]["VolumeMelt(m³)"] = compo_melt_all_volume
        t["Melt"]["NameUnitMelting"] = compo_melt_all_unit
    else
        t["Melt"]["CompositionMelt(wt%)"] = convert(Array{Float32}, zeros(11))
        t["Melt"]["IndexGridMelt"] = [0]
        t["Melt"]["TimeMelt(Myr)"] = [0]
        t["Melt"]["VolumeMelt(m³)"] = [0]
        t["Melt"]["NameUnitMelting"] = [""]
    end
    t["Rock"]["AssemblageName"] = magemin_interact.name_assemblage
    t["Rock"]["AssemblageMode(mol)"] = convert(Array{Float32}, magemin_interact.mode_assemblage)
    t["Rock"]["AssemblageMode(vol)"] = convert(Array{Float32}, magemin_interact.mode_assemblage_vol)
end


"""
    hdf5_initial_conditions(integrator, path_hdf5)

This function saves the initial conditions of the model to an HDF5 file. It takes the `integrator` object and the path to the HDF5 file `path_hdf5` as arguments. It creates groups and datasets in the HDF5 file and sets their attributes.
"""
function hdf5_initial_conditions(integrator, path_hdf5)

    @unpack grid, thermal_param, compositions, physical_prop, magemin_interact = integrator.p
    @unpack depthc, depthc_ini, Δz, nz, tfinal, Lz = grid
    @unpack T0, α, Cp, k, A, ρ, Ra, Nu = thermal_param
    @unpack Plith = physical_prop
    @unpack compo_solid, compo_solid_matrix, compo_solid_name = compositions
    @unpack compo_rock_wt, compo_i_g, total_volume_melt, compo_melt_wt_all, compo_melt_wt_all_index, compo_melt_all_time, oxides_name, crystal_fractionation, frac_crystallisation, depth_crystal_fractionation = magemin_interact

    compo_solid_matrix .= hcat(compo_solid...)'

    convert_compo_solid_to_wt_pct!(compo_rock_wt, compo_solid_matrix, compo_i_g)


    h5open(path_hdf5, "w") do file
        g = create_group(file, "MeltingCrust") # create a group

        # Set attributes for group g
        attrs_values = (
            ("Lz(m)", Lz),
            ("Dz(m)", Δz),
            ("Nz", nz),
            ("TotalTime(Myr)", tfinal / (365.25 * 24 * 3600 * 1e6)),
            ("Depth(m)", collect(depthc)),
            ("OxidesNames", oxides_name),
            ("FractionalCrystallisation", crystal_fractionation),
            ("FracCrystallization", frac_crystallisation),
            ("DepthFractionCrystallisation", depth_crystal_fractionation)
        )

        for (attr, value) in attrs_values
            attributes(g)[attr] = value
        end

        t = create_group(file, "MeltingCrust/t$(lpad("0", 6, "0"))") # create a group

        # Create groups and set attributes for each group
        groups = ("Density(kg.m⁻³)", "Cp(J.kg⁻¹.K⁻¹)", "k(W.m⁻¹.K⁻¹)", "α(m².s⁻¹)")
        for group in groups
            grp = create_group(t, group)
            attributes(grp)["DataType"] = "Scalar"
            attributes(grp)["Center"] = "Node"
        end

        # Create datasets and set attributes for each dataset
        t["Density(kg.m⁻³)"]["Density(kg.m⁻³)"] = ρ
        t["Cp(J.kg⁻¹.K⁻¹)"]["Cp(J.kg⁻¹.K⁻¹)"] = convert(Array{Float32}, Cp)
        t["k(W.m⁻¹.K⁻¹)"]["k(W.m⁻¹.K⁻¹)"] = convert(Array{Float32}, k)
        t["α(m².s⁻¹)"]["α(m².s⁻¹)"] = convert(Array{Float32}, α)

        set_attributes_and_groups(t, 0, 0, magemin_interact, integrator, Plith, compo_solid_matrix, compo_solid_name)
    end

end

"""
    hdf5_timestep(integrator, Δt, tcurrent, path_hdf5)

This function saves the state of the model at a specific timestep to an HDF5 file. It takes the `integrator` object, time step `Δt`, current time `tcurrent`, and the path to the HDF5 file `path_hdf5` as arguments. It creates groups and datasets in the HDF5 file and sets their attributes.
"""
function hdf5_timestep(integrator, Δt, tcurrent, path_hdf5)

    @unpack grid, thermal_param, compositions, physical_prop, magemin_interact = integrator.p
    @unpack depth, Δz, nz, tfinal, Lz = grid
    @unpack T0, α, Cp, k, A, ρ, Ra, Nu = thermal_param
    @unpack Plith = physical_prop
    @unpack compo_solid, compo_solid_matrix, compo_solid_name = compositions
    @unpack compo_rock_wt, compo_i_g, total_volume_melt, compo_melt_wt_all, compo_melt_wt_all_index, compo_melt_all_time = magemin_interact

    compo_solid_matrix .= hcat(compo_solid...)'

    # convert compo_rock_wt from mol% to wt%
    convert_compo_solid_to_wt_pct!(compo_rock_wt, compo_solid_matrix, compo_i_g)

    h5open(path_hdf5, "r+") do file
        # output the number of group in the HDF5
        n = length((file["MeltingCrust"]))

        t = create_group(file, "MeltingCrust/t$(lpad(string(n), 6, "0"))") # create a group

        set_attributes_and_groups(t, tcurrent, Δt, magemin_interact, integrator, Plith, compo_solid_matrix, compo_solid_name)
    end
end


"""
    save_data(integrator)

Callback function used to save data produced by the model to an HDF5 file at a specific timestep.
This function checks the current time of the `integrator` and calls the appropriate function to save the data.
"""
function save_data(integrator)

    @unpack path_save = integrator.p

    if integrator.t ≠ 0.0
        hdf5_timestep(integrator, integrator.dt, integrator.t, path_save)
        # @info "Data saved at $(round((integrator.t / (365.25 * 24 * 3600 * 1e6)), digits=2)) Myr."
    elseif integrator.t == 0.0
        hdf5_initial_conditions(integrator, path_save)
        # @info "Data saved at $(round((integrator.t / (365.25 * 24 * 3600 * 1e6)), digits=2)) Myr."
    end

    # @unpack geoparam, grid, thermal_param, physical_prop, magemin_interact = integrator.p
    # @unpack z, zc = grid
    # @unpack k, Cp, ρ, Ql, A, k_conv, q_conv, Nu, convection, advection, input_gabbro, heat_flux_moho, neumann_bottom = thermal_param
    # @unpack ϕ, dϕdT = physical_prop
    # MatParam, phases = geoparam
    # @unpack depth_min_thermo, depth_max_thermo= magemin_interact


    # PT_points_T = [799.4632313472895,809.9469940955452,815.1888754696728,878.0914519592058,880.7123926462697,757.5281803542675,707.7303073000539,600.2717391304349,583.2356246645198]
    # PT_points_depth = [17.242374103771585,18.87888586610666,20.26681356327691,21.571879905392223,23.47769107165585,16.84878266726062,15.502285647617837,12.581422881931188,11.607801960035639]
    # P_uncertainties = [0.6,0.6,0.6,0.6,0.6,0.6,2.5,2.5,2.5,]
    # T_uncertainties = [55,55,55,55,55,55,55,55,55]

    # layout = @layout [a b]
    # a = plot(grid.depthc.*1e-3, integrator.u .- 273.15, label = "Initial conditions", dpi=400, title = "Total time: $(round((integrator.t ./ (1e6*365.25*24*3600)); digits=2)) Ma")
    # scatter!(a, PT_points_depth, PT_points_T, yerror=T_uncertainties, xerror=P_uncertainties, label = "PT points", marker = :circle, markersize = 3, legend = :topleft)
    # b = plot(grid.depthc.*1e-3, Cp, label = "Cp", dpi=400, title = "Cp", ylims=(0, 10000))

    # display(plot(a,b, layout=layout))

end
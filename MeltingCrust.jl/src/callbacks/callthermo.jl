
function clear_array(array::Array{T}) where T
    if !isempty(array)
        for i in axes(array, 1)
            if !isempty(array[i])
                array[i] = eltype(array)[]
            end
        end
    end
end

function clear_array(array::Array{T}) where T <: Number
    if !iszero(array)
        for I in CartesianIndices(array)
            if !iszero(array[I])
                array[I] = 0
            end
        end
    end
end

function clear_array(array::Array{Bool})
    if !isempty(array)
        for i in axes(array, 1)
            if !isempty(array[i])
                array[i] = false
            end
        end
    end
end

function minimization_calls_magemin!(integrator)

    @unpack grid, compositions, physical_prop, thermal_param, magemin_interact, data_MAGEMin = integrator.p
    @unpack depth, depthc, depthc_ini, nz, z, zc = grid
    @unpack depth_min_thermo, depth_max_thermo, miniz_calls, ΔT, time_melt, time_volume, oxides_name, vol_melt_extract, vol_frac_system_shrinking_tot, compo_melt_wt, compo_melt_wt_all, compo_melt_wt_all_index, compo_melt_all_time, compo_melt_wt_dry_i, compo_i_g, compo_melt_i_wt, flux_melting, melt_prev, compo_melt_all_volume, compo_melt_all_unit, ms_melt_check, grt_frac, grt_frac_prev, grt_frac_remove, crd_frac, crd_frac_prev, crd_frac_remove = magemin_interact
    @unpack compo_solid, compo_solid_name = compositions
    @unpack Plith = physical_prop
    @unpack ρ, Cp = thermal_param


    clear_array(magemin_interact.updated)
    clear_array(magemin_interact.compo_melt_wt)

    update_temperature_and_minimization_calls!(integrator, magemin_interact)

    # check if magemin_interact.updated is all false
    if !iszero(magemin_interact.updated)

        if flux_melting
            water_flux_melting!(compo_solid, magemin_interact)
        end

        out = perform_minimization(magemin_interact, Plith, integrator, compo_solid, data_MAGEMin, oxides_name)

        for (i, index) = enumerate(findall(magemin_interact.updated))
            @unpack frac_M, bulk, bulk_M, frac_S, bulk_S, ph, ph_frac, ph_frac_vol, ph_type, ph_id, SS_vec, PP_vec, oxides =  out[i];

            # ρ[index] = out[i].rho  # update density for thermal model
            # Cp[index] = out[i].s_cp[1]  # update heat capacity for thermal model


            find_ox_index!(magemin_interact.el_index_MAGEMin, magemin_interact.el_index_meltingcrust, magemin_interact.oxides_name, oxides)

            # remove H2O if presents
            extract_H2O_MAGEMin!(bulk, out[i])
            # remove melt and compute dry melt composition


            # normalise grt_frac without liq if present
            if any(occursin.("g"), out[i].ph)
                if any(occursin.("liq", out[i].ph))
                    id_grt = findfirst(out[i].ph .== "g")
                    id_liq = findfirst(out[i].ph .== "liq")
                    grt_frac[index] = ph_frac[id_grt] / (1 - ph_frac[id_liq])
                # else
                #     id_grt = findfirst(out[i].ph .== "g")
                #     grt_frac[index] = ph_frac[id_grt]
                end
            end

            # normalise crd_frac without liq if present
            if any(occursin.("cd"), out[i].ph)
                if any(occursin.("liq", out[i].ph))
                    id_crd = findfirst(out[i].ph .== "cd")
                    id_liq = findfirst(out[i].ph .== "liq")
                    crd_frac[index] = ph_frac[id_crd] / (1 - ph_frac[id_liq])
                # else
                #     id_crd = findfirst(out[i].ph .== "cd")
                #     crd_frac[index] = ph_frac[id_crd]
                end
            end

            # deal with melt if higher than 7 vol%
            if frac_M >= 0.07  && ms_melt_check[index] == true

                vol_melt_extract[index], vol_frac_system_shrinking_tot[index], vol_frac_system_shrinking = extract_melt_MAGEMin!(bulk, vol_frac_system_shrinking_tot[index], magemin_interact, out[i], grid, index)
                # store composition of melt with water in wt%

                # remove volume corresponding to the melt extracted
                # reduction_Δz = (depth[index] - depth[index+1]) * (vol_frac_system_shrinking)  # length of the rock extracted
                # depth[index:-1:1] .-= reduction_Δz  # remove this length on the depth bellow
                # depthc[index:-1:1] = (depth[index:-1:1] + depth[index+1:-1:2]) / 2  # recenter the center of the cell
                # # depthc[index-1:-1:1] .+= reduction_Δz / 2 # recenter the center of the bellow cells
                # @show vol_frac_system_shrinking_tot[index], reduction_Δz
                # @show depth[index+1:-1:index-2]
                # @show depthc[index+1:-1:index-2]
                # @show index, Plith[index], integrator.u[index].-273.15

                #! call  crystal fractionation here
                if magemin_interact.crystal_fractionation
                    ph_frac_tot = process_crystal_fractionation!(bulk_M, depthc_ini, Plith, integrator, magemin_interact, data_MAGEMin, oxides_name)
                    vol_melt_extract[index] *= (1 - ph_frac_tot)  # remove the fraction of crystal from melt
                end

                # add peritectic garnet to melt
                if magemin_interact.grt_frac_remove > 0
                    # add the composition of the grt to bulk_M
                    if any(occursin.("g"), out[i].ph)
                        # if grt_frac[index] > grt_frac_prev[index] # if grt is increasing
                        id_grt = findfirst(out[i].ph .== "g")
                        @show grt_frac[index]
                        @show bulk_M
                        bulk_M .= bulk_M .+ grt_frac_remove .* grt_frac[index] .* out[i].SS_vec[id_grt].Comp
                        @show bulk_M

                        # remove the composition of the grt from bulk
                        bulk .= bulk .- grt_frac_remove .* grt_frac[index] .* out[i].SS_vec[id_grt].Comp

                        # end
                    end
                end

                # add peritectic cordierite to melt
                if magemin_interact.crd_frac_remove > 0
                    # add the composition of the grt to bulk_M
                    if any(occursin.("cd"), out[i].ph)
                        # if grt_frac[index] > grt_frac_prev[index] # if grt is increasing
                        id_crd = findfirst(out[i].ph .== "cd")
                        @show bulk_M
                        bulk_M .= bulk_M .+ magemin_interact.crd_frac_remove .* out[i].SS_vec[id_crd].Comp
                        @show bulk_M

                        # remove the composition of the grt from bulk
                        bulk .= bulk .- crd_frac_remove .* crd_frac[index] .* out[i].SS_vec[id_crd].Comp

                        # end
                    end
                end


                magemin_interact.total_volume_melt += vol_melt_extract[index]

                bulk_M_ind = @view bulk_M[magemin_interact.el_index_meltingcrust]

                # normalize bulk_M from mol% to wt%
                convert_compo_to_wt_pct!(compo_melt_i_wt, bulk_M_ind, compo_i_g)

                compo_melt_wt[index,:] .= compo_melt_i_wt ./ sum(compo_melt_i_wt)

                # save melt composition in wt%
                compo_melt_wt_index = index

                # for j in compo_melt_wt_index
                #     # convert wt% to dry wt%
                #     convert_compo_to_wt_dry!(compo_melt_wt_dry_i, compo_melt_wt[j,:])

                #     push!(magemin_interact.compo_melt_wt_all, compo_melt_wt[j,:])
                #     push!(magemin_interact.compo_melt_wt_dry_all, compo_melt_wt_dry_i[:])
                # end

                # convert wt% to dry wt%
                convert_compo_to_wt_dry!(compo_melt_wt_dry_i, compo_melt_wt[index,:])

                push!(magemin_interact.compo_melt_wt_all, compo_melt_wt[index,:])
                push!(magemin_interact.compo_melt_wt_dry_all, compo_melt_wt_dry_i[:])

                append!(compo_melt_wt_all_index, compo_melt_wt_index)
                # time_melt = Vector{Float64}(undef,length(compo_melt_wt_index))
                # vol = Vector{Float64}(undef,length(compo_melt_wt_index))

                # time_melt .= integrator.t
                append!(compo_melt_all_time, integrator.t)
                append!(compo_melt_all_volume, vol_melt_extract[index])
                push!(compo_melt_all_unit, compo_solid_name[index])

                # update melt_prev at this index to 1
                melt_prev[index] = 1

                if magemin_interact.ms_melt_only
                    # check if out.ph contains "mu" and "liq"
                    if any(occursin.("mu", out[i].ph)) && any(occursin.("liq", out[i].ph))
                        ms_melt_check[index] = true
                    else
                        ms_melt_check[index] = false
                    end
                end
            else
                # update melt_prev at this index to 0
                melt_prev[index] = 0
            end

            # update new compo from bulk
            compo_solid[index][magemin_interact.el_index_MAGEMin] .= bulk

            # save part
            ph_length = axes(ph,1)


            # if there is duplicate in ph, then print it
            # append the name of the phases

            magemin_interact.name_assemblage[index,:] .= fill("", 11)
            magemin_interact.mode_assemblage[index,:] .= fill(0.0, 11)
            magemin_interact.mode_assemblage_vol[index,:] .= fill(0.0, 11)

            magemin_interact.name_assemblage[index,ph_length] .= ph


            magemin_interact.mode_assemblage[index,ph_length] .= ph_frac
            magemin_interact.mode_assemblage_vol[index,ph_length] .= ph_frac_vol

            # reset composition assemblage matrix if magemin was called
            magemin_interact.compo_assemblage[magemin_interact.el_index_MAGEMin,(index-1)*11:((index-1)*11+11)] .= 0

            for j in axes(ph,1)
                if j ≤ length(SS_vec)
                    magemin_interact.compo_assemblage[magemin_interact.el_index_MAGEMin,(index-1)*11 + j] .= SS_vec[j].Comp
                else
                    magemin_interact.compo_assemblage[magemin_interact.el_index_MAGEMin, (index-1)*11 + j] .= PP_vec[j-length(SS_vec)].Comp
                end
            end

            # save grt_frac in grt_frac_prev
            # if any(occursin.("g"), out[i].ph)
            #     grt_frac_prev[index] = grt_frac[index]
            #     grt_frac[index] = 0
            # end

        end
    end

    # renormalize bulk composition
    for i in axes(compo_solid,1)
        compo_solid[i,:] ./= sum(compo_solid[i])
    end


    l = @layout [a b c]

    p1 = plot(ΔT, depthc.*0.001,  xlabel="ΔT", ylabel= "Depth (km)", label="", yflip = true)
    p2 = plot(miniz_calls, depthc.*0.001,  xlabel="Total number of \nminimization calls", ylabel= "Depth (km)", label="", yflip = true)
    p3 = plot(integrator.u .- 273.15, depthc.*0.001,  xlabel="Temperature (°C)", ylabel= "Depth (km)", yflip = true, label="")


    display(plot(p1, p2, p3, layout=l, plot_title="Total time = $(integrator.t / (365.25 * 24 * 3600 * 1e6)) Myr", size= (800, 800)))

end

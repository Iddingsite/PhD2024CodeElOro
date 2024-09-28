using Parameters

function find_ox_index!(el_index_MAGEMin, el_index_meltingcrust, oxides_name, oxides_MAGEMin)

    # check if they are already in the good order from previous calls
    for i in axes(oxides_name,1)
        el_index_MAGEMin[i] = findfirst(==(oxides_MAGEMin[i]), oxides_name)
    end

    # do the same but for the opposite
    for i in axes(oxides_MAGEMin,1)
        el_index_meltingcrust[i] = findfirst(==(oxides_name[i]), oxides_MAGEMin)
    end

end

function extract_H2O_MAGEMin!(bulk, out)

    # check if ph contains "H2O"
    if any(occursin.("H2O", out.ph))
        id_H2O = findfirst(out.ph .== "H2O")
        bulk .= bulk .- out.bulk_F .* out.ph_frac[id_H2O]
        # put frac of H2O to 0
        out.ph_frac[id_H2O] = 0
        out.ph_frac_vol[id_H2O] = 0
        # renormalise to 100
        out.ph_frac .= out.ph_frac ./ sum(out.ph_frac)
        out.ph_frac_vol .= out.ph_frac_vol ./ sum(out.ph_frac_vol)
    end

end

function extract_melt_MAGEMin!(bulk, vol_frac_system_shrinking_tot, magemin_interact, out, grid, index)
    @unpack frac_M, bulk_M, bulk_M_wt, bulk_wt =  out;
    @unpack Δz, Lx = grid
    @unpack grt_frac, grt_frac_prev, grt_frac_remove = magemin_interact


    # extract only if melt is more than 7 vol percent
    vol_remove = 1 - 1/(100*frac_M)  # keep 1% of melt
    if sum(bulk_M) != 1
        bulk_M = bulk_M ./ sum(bulk_M)
    end  # normalize to 1

    # slope = 0

    # grt_present = any(occursin.("g", out.ph))
    # # compute slope of grt_frac and grt_frac_prev when there is melt
    # if grt_present
    #     if grt_frac[index] > grt_frac_prev[index] # if grt is increasing
    #         slope = 1
    #     else
    #         slope = -1
    #     end
    # end

    bulk .= bulk .- vol_remove .* bulk_M .* frac_M  # remove melt from bulk

    # # remove grt from bulk
    # if slope == 1
    #     # get composition of grt in out.SS_vec
    #     id_grt = findfirst(out.ph .== "g")
    #     bulk .= bulk .- grt_frac_remove .* out.SS_vec[id_grt].Comp .* grt_frac[index]
    # end


    # find index of liq
    id_liq = findfirst(out.ph .== "liq")

    # store volumed extracted
    # vol_l = molar_vol_l * mole_l * volume_tot / molar_vol_tot / mole_tot
    vol_tot = Lx * Δz * (1-vol_frac_system_shrinking_tot)  # volume of the system
    vol_frac_system_shrinking = out.ph_frac_vol[id_liq] * vol_remove  # volume fraction extracted while keeping 1% of melt
    vol_melt_extract_i = vol_frac_system_shrinking * vol_tot  # volume of melt extracted

    vol_frac_system_shrinking_tot = vol_melt_extract_i / (vol_tot) + vol_frac_system_shrinking_tot  # vol extract / total vol + percent of volume already extracted -> total fraction of thinning or total fraction of melt extracted

    return vol_melt_extract_i, vol_frac_system_shrinking_tot, vol_frac_system_shrinking
end

function update_temperature_and_minimization_calls!(integrator, magemin_interact)
    @unpack ΔT, depth_min_thermo, depth_max_thermo, miniz_calls = magemin_interact
    @unpack depthc_ini = integrator.p.grid

    for i in axes(integrator.u,1)
        ΔT_timestep = integrator.u[i] - integrator.uprev[i]
        if ΔT_timestep > 0
            ΔT[i] += ΔT_timestep
        end

        if ΔT[i] > 5
            if depthc_ini[i] > depth_min_thermo && depthc_ini[i] < depth_max_thermo
                magemin_interact.updated[i] = true
                miniz_calls[i] += 1
                ΔT[i] = 0
            end
        end
    end
end


function process_crystal_fractionation!(bulk_M, depthc_ini, Plith, integrator, magemin_interact, data_MAGEMin, oxides_name)
    @unpack depth_crystal_fractionation, frac_crystallisation = magemin_interact

    bulk_M_ind = @view bulk_M[magemin_interact.el_index_meltingcrust]

    # find index of depth_crystal_fractionation in depth
    depth_index = findmin(abs.(depthc_ini .- depth_crystal_fractionation))[2]

    Plith_frac = round(Plith[depth_index] * 1e-8, digits=2)
    T_frac = round(integrator.u[depth_index] - 273.15, digits=2)

    # bulk_M_frac = zeros(size(bulk_M))
    # bulk_M_frac[magemin_interact.el_index_MAGEMin] .= bulk_M
    bulk_M_ind .= max.(bulk_M_ind, 0.0001)
    bulk_M_ind .= round.(bulk_M_ind, digits=4)
    # @show bulk_M, Plith_frac, T_frac
    output_melt = single_point_minimization(Plith_frac, T_frac, data_MAGEMin, X=bulk_M_ind, Xoxides=oxides_name, sys_in="mol")

    # find all position of "fsp" in output_melt.ph
    fsp_index = findall(x -> x == "fsp", output_melt.ph)
    id_qz = findfirst(output_melt.ph .== "q")
    id_PP_qz = output_melt.ph_id[id_qz] + 1  # +1 because index starts at 0 (C)

    ph_frac_tot = 0  # fraction extracted from the mineral assemblage

    # if there is fsp in the phase, extract composition from output_melt.SS_vec
    if !isempty(fsp_index)
        for j in fsp_index
            # get index for K2O, Na2O and CaO in output_melt.oxides
            k2o_index = findfirst(output_melt.oxides .== "K2O")
            na2o_index = findfirst(output_melt.oxides .== "Na2O")
            cao_index = findfirst(output_melt.oxides .== "CaO")

            # extract composition of fsp
            fsp_comp = output_melt.SS_vec[j].Comp

            # convert k2o, na2o and cao to cation of K, Na and Ca
            K = fsp_comp[k2o_index] * 2
            Na = fsp_comp[na2o_index] * 2
            Ca = fsp_comp[cao_index]

            # calculate K / (K + Na + Ca)
            kfs = K / (K + Na + Ca)

            # get fraction of fsp in output_melt.ph_frac
            fsp_frac = output_melt.ph_frac[j]
            fsp_frac_vol = output_melt.ph_frac_vol[j]

            # remove fsp from bulk and save the vol fraction of fsp
            if kfs < 0.10  # if kfs < 10%, then it is Plagioclase
                output_melt.bulk .= output_melt.bulk .- fsp_frac .* fsp_comp .* frac_crystallisation[1]

                ph_frac_tot += fsp_frac_vol .* frac_crystallisation[1]
            else
                output_melt.bulk .= output_melt.bulk .- fsp_frac .* fsp_comp .* frac_crystallisation[2]

                ph_frac_tot += fsp_frac_vol .* frac_crystallisation[2]
            end

            # replace values bellow 0 by 0
            output_melt.bulk .= max.(output_melt.bulk, 0.0001)
            # renormalize to 100
            output_melt.bulk .= output_melt.bulk ./ sum(output_melt.bulk)
        end
    end

    if !isempty(id_PP_qz)
        qz_comp = output_melt.PP_vec[id_PP_qz].Comp
        qz_frac = output_melt.ph_frac[id_qz]
        qz_frac_vol = output_melt.ph_frac_vol[id_qz]

        output_melt.bulk .= output_melt.bulk .- qz_frac .* qz_comp .* frac_crystallisation[3]
        ph_frac_tot += qz_frac_vol .* frac_crystallisation[3]

        # replace values bellow 0 by 0
        output_melt.bulk .= max.(output_melt.bulk, 0.0001)

        # renormalize to 100
        output_melt.bulk .= output_melt.bulk ./ sum(output_melt.bulk)
    end

    bulk_M .= output_melt.bulk

    bulk_M .= round.(bulk_M, digits=4)

    return ph_frac_tot
end

function set_value_and_normalize!(array, index, new_value)
    # Set the desired value
    array[index] = new_value

    # Calculate the sum of the other values
    other_sum = sum(array) - new_value

    # Calculate the factor to multiply the other values by
    factor = (1 - new_value) / other_sum

    # Normalize the other values
    for i in axes(array,1)
        if i != index
            array[i] *= factor
        end
    end
end

function water_flux_melting!(compo_solid, magemin_interact)
    for i in axes(magemin_interact.updated, 1)
        if magemin_interact.updated[i] && magemin_interact.melt_prev[i]
            # renormalize compo_solid to 1
            compo_solid[i][:] ./= sum(compo_solid[i])
            # update H2O in compo_solid to the value of h2o_keep if it is lower than h2o_keep
            # if compo_solid[i][end-1] < magemin_interact.h2o_keep
                # renormalise compo_solid to 1 by forcing H2O to h2o_keep
                index = length(compo_solid[i])-1
                compo_solid_i = @view compo_solid[i][:]
                water_to_add = magemin_interact.h2o_keep + compo_solid_i[end-1]
                set_value_and_normalize!(compo_solid_i, index, water_to_add)

                # renormalise the rest to be equal to 1-h2o_keep
                # from index 1:end-2
            # end
        end
    end
end


function perform_minimization(magemin_interact, Plith, integrator, compo_solid, data_MAGEMin, oxides_name)

    P_i = round.(Plith[magemin_interact.updated].*1e-8; digits=2)
    T_i = round.(integrator.u[magemin_interact.updated].-273.15; digits=2)


    # count_up = count(magemin_interact.updated)
    index_update = findall(magemin_interact.updated)

    # if count_up > 1
    #     compo_solid_i = @view compo_solid[index_update]
    # else # count_up is 1
    #     index_updated = index_update[1]
    #     # prevent view of dimension 0
    #     compo_solid_i = @view compo_solid[index_updated:index_updated][:]
    # end

    # @show compo_solid[index_update]
    # compo_solid_i = copy([round.(sub_array; digits=4) for sub_array in compo_solid[index_update]])
    # @show compo_solid_i, T_i, P_i

    compo_solid_i = [(max.(sub_array, 0.0001)) for sub_array in compo_solid[index_update]]

    # @show compo_solid_i, T_i, P_i
    # println("MAGEMin: start minimization")
    output=Vector{MAGEMin_C.gmin_struct{Float64, Int64}}(undef,length(index_update))
    output .= multi_point_minimization(P_i, T_i, data_MAGEMin, X=compo_solid_i, Xoxides=oxides_name, sys_in="mol", scp=0)
    # println("MAGEMin: finish minimization")

    return output
end

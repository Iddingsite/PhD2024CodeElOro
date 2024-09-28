function initial_Plith!(Plith, ρ, z)
    # Initial Plith
    for i in axes(Plith,1)
        if i == 1
            Plith[i] = 0
        else
            Plith[i] = Plith[i-1] + ρ[i] * 9.81 * (z[i] - z[i-1])
        end
    end

    Plith .= reverse(Plith)
end

function initial_T!(T0, zc)
    for i in axes(zc,1)
        T0[i] = 0.025*zc[i] + 0 + 273.15
    end
    reverse!(T0, dims=1)
end

function initial_Nu!(Nu, Ra, zc)
    # Initial Nu
    for i in axes(zc,1)
        if reverse(zc)[i] >= 17_000 && reverse(zc)[i] <= 25_000
            Nu[i] = 0.120*Ra^(1/3)
        end
    end
end

function layers_compo!(compo_solid, compo_solid_name, compo_pelite, compo_arkose, compo_quartzite)
    counter_compo = 1
    for i in axes(compo_solid,1)
        if counter_compo == 1
            compo_solid_name[i] = "pelite"
        elseif counter_compo == 2
            compo_solid_name[i] = "semi-pelite"
        elseif counter_compo == 3
            compo_solid_name[i] = "quartzite"
        end
        for k in axes(compo_solid[1],1)
            if counter_compo == 1
                compo_solid[i][k] = compo_pelite[k]
            elseif counter_compo == 2
                compo_solid[i][k] = compo_arkose[k]
            elseif counter_compo == 3
                compo_solid[i][k] = compo_quartzite[k]
            end
        end
        counter_compo+=1
        if counter_compo == 4
            counter_compo = 1
        end
    end
end


function initial_call_MAGEMin!(compo_solid, grid_model, physical_prop, thermal_param, magemin_interact, T0, data_MAGEMin)
    @unpack depth = grid_model
    @unpack depth_min_thermo, depth_max_thermo, miniz_calls, ΔT, oxides_name = magemin_interact
    @unpack Plith = physical_prop
    @unpack ρ, Cp = thermal_param

    # MeltingCrust.clear_array(magemin_interact.name_assemblage)
    # MeltingCrust.clear_array(magemin_interact.mode_assemblage)
    MeltingCrust.clear_array(magemin_interact.updated)


    for i in axes(T0,1)
        if depth[i] > depth_min_thermo && depth[i] < depth_max_thermo
            magemin_interact.updated[i] = true
        end
    end

    if length(magemin_interact.updated) > 1
        compo_solid_i = @view compo_solid[magemin_interact.updated][:]
    elseif length(magemin_interact.updated) == 0
        compo_solid_i = @view compo_solid[magemin_interact.updated:magemin_interact.updated][:]
    end

    P_i = Plith[magemin_interact.updated]
    T_i = T0[magemin_interact.updated]

    # out = multi_point_minimization(P_i.*1e-8, T_i.-273.15, data_MAGEMin, X=compo_solid_i, Xoxides=oxides_name, sys_in="mol", scp=1);
    out = multi_point_minimization(P_i.*1e-8, T_i.-273.15, data_MAGEMin, X=compo_solid_i, Xoxides=oxides_name, sys_in="mol", scp=0);

    for (i, index) = enumerate(findall(magemin_interact.updated))

        @unpack frac_M, bulk, bulk_M, frac_S, bulk_S, ph, ph_frac, ph_type, ph_id, SS_vec, PP_vec, oxides =  out[i];

        # ρ[index] = out[i].rho

        # Cp[index] = out[i].s_cp[1]

        find_ox_index!(magemin_interact.el_index_MAGEMin, magemin_interact.el_index_meltingcrust, magemin_interact.oxides_name, oxides)

        # remove H2O if presents
        extract_H2O_MAGEMin!(bulk, out[i])
        # append the name of the phases
        # remove melt and compute dry melt composition
        # deal with melt if higher than 7 vol%
        if frac_M >= 0.07

            vol_melt_extract[index] = extract_melt_MAGEMin!(bulk, out[i], grid)
            # store composition of melt with water in wt%

            magemin_interact.total_volume_melt += vol_melt_extract[index]

            compo_melt_wt[index,magemin_interact.el_index_MAGEMin] .= bulk_M ./ sum(bulk_M)


            # save melt composition in wt%
            compo_melt_wt_index = findall(!iszero, compo_melt_wt[:,1])

            for i in compo_melt_wt_index
                # convert wt% to dry wt%
                convert_compo_to_wt_dry!(compo_melt_wt_dry_i, compo_melt_wt[i,:])

                push!(magemin_interact.compo_melt_wt_all, compo_melt_wt[i,:])
                push!(magemin_interact.compo_melt_wt_dry_all, compo_melt_wt_dry_i[:])
            end

            append!(compo_melt_wt_all_index, compo_melt_wt_index)
            time_melt = Vector{Float64}(undef,length(compo_melt_wt_index))
            time_melt .= integrator.t
            append!(compo_melt_all_time, time_melt)
        end

        # update new compo from bulk
        compo_solid[index][magemin_interact.el_index_MAGEMin] .= bulk

        ph_length = axes(ph,1)

        magemin_interact.name_assemblage[index,:] .= fill("", 11)
        magemin_interact.mode_assemblage[index,:] .= fill(0.0, 11)

        magemin_interact.name_assemblage[index,ph_length] .= ph
        magemin_interact.mode_assemblage[index,ph_length] .= ph_frac

        # reset composition assemblage matrix if magemin was called
        magemin_interact.compo_assemblage[magemin_interact.el_index_MAGEMin,(index-1)*11:((index-1)*11+11)] .= 0

        for j in axes(ph,1)
            if j ≤ length(SS_vec)
                magemin_interact.compo_assemblage[magemin_interact.el_index_MAGEMin,(index-1)*11 + j] .= SS_vec[j].Comp
            else
                magemin_interact.compo_assemblage[magemin_interact.el_index_MAGEMin, (index-1)*11 + j] .= PP_vec[j-length(SS_vec)].Comp
            end
        end


    end


end
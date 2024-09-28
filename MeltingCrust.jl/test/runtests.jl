using Test
using MeltingCrust

@testset "MAGEMin wrapper" begin

    data = Initialize_MAGEMin("mp", verbose=true);
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "MnO"; "H2O"; "O"]
    bulk_in_forshaw = [70.48; 25.43/2; 0.77; 3.95; 6.3; 5.54/2; 2.94/2; 0.75; 0.07/2; 100/2; 0]

    P = 4.0  # in kbar
    T = 500.0  # in °C

    sys_in = "mol"

    out     = single_point_minimization(P, T, data, X=bulk_in_forshaw, Xoxides=Xoxides)

    @test any(out.ph .== "H2O") == true  # check if H2O is present
    @test any(out.ph .== "liq") == false  # check if melt is present
end

@testset "MAGEMin extraction fluids" begin

    data = Initialize_MAGEMin("mp", verbose=false);
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "MnO"; "H2O"; "O"]
    bulk_in_forshaw = [70.48; 25.43/2; 0.77; 3.95; 6.3; 5.54/2; 2.94/2; 0.75/2; 0.07; 100/2; 0.1]

    P = 4.0
    T = 500.0

    sys_in = "mol"

    out     = single_point_minimization(P, T, data, X=bulk_in_forshaw, Xoxides=Xoxides)
    @unpack frac_M, bulk, bulk_M, frac_S, bulk_S, ph, ph_frac, ph_type, ph_id, SS_vec, PP_vec =  out;

    id_H2O = findfirst(ph .== "H2O")

    new_bulk = bulk .- out.bulk_F .* ph_frac[id_H2O]

    # put frac of H2O to 0
    out.ph_frac[id_H2O] = 0
    # renormalise to 100
    out.ph_frac .= out.ph_frac ./ sum(out.ph_frac)

    # test function now
    out2 = single_point_minimization(P, T, data, X=bulk_in_forshaw, Xoxides=Xoxides)

    extract_H2O_MAGEMin!(out2.bulk, out2)

    @test new_bulk[end] ≈ 0.08612 rtol=1e-4
    @test out2.bulk[end] ≈ 0.08612 rtol=1e-4
end


@testset "MAGEMin extraction of melt" begin


    data = Initialize_MAGEMin("mp", verbose=false);
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "MnO"; "H2O"; "O"]
    bulk_in_forshaw = [70.48; 25.43/2; 0.77; 3.95; 6.3; 5.54/2; 2.94/2; 0.75/2; 0.07; 100/2; 0.1]

    P = 8.0
    T = 700.0

    sys_in = "mol"

    out     = single_point_minimization(P, T, data, X=bulk_in_forshaw, Xoxides=Xoxides)

    grid_model = CreateGrid(nz=100, Lz=28u"km", tfinal=10u"Myr")

    magemin_interact = MAGEMinInteraction(grid=grid_model, crystal_fractionation=false, frac_crystallisation=[0, 0, 0.0], flux_melting=false, h2o_keep= 0.0, ms_melt_only=false)

    @unpack frac_M, bulk, bulk_M, frac_S, bulk_S, ph, ph_frac, ph_type, ph_id, SS_vec, PP_vec, bulk_wt =  out;

    bulk_test = copy(bulk)

    vol_frac_system_shrinking_tot = 0

    # extract only if melt is more than 7 vol percent
    if frac_M >= 0.07
        @unpack Δz, Lx = grid_model
        vol_remove = (100-1/frac_M) / 100

        bulk_test .= bulk .- vol_remove .* bulk_M .* frac_M

        # store volumed extracted for a 2D surface
        #! check that it works within the function
            # find index of liq
        id_liq = findfirst(x -> x == "liq", out.ph)

        # store volumed extracted
        # vol_l = molar_vol_l * mole_l * volume_tot / molar_vol_tot / mole_tot
        vol_tot = Lx * Δz * (1-vol_frac_system_shrinking_tot)  # volume of the system
        vol_frac_system_shrinking = out.ph_frac_vol[id_liq]  * vol_remove # volume fraction extracted
        vol_melt_extract_i = vol_frac_system_shrinking * vol_tot  # volume of melt extracted

        vol_frac_system_shrinking_tot = vol_melt_extract_i / (vol_tot) + vol_frac_system_shrinking_tot  # vol extract / total vol + percent of volume already extracted -> total fraction of thinning or total fraction of melt extracted

    end

    @test sum(bulk_test) ≈ 0.7815 rtol=1e-4
    @test vol_remove ≈ 0.9562  rtol=1e-4
    @test vol_melt_extract_i ≈ 65.0696 rtol=1e-4

    out2     = single_point_minimization(P, T, data, X=bulk_in_forshaw, Xoxides=Xoxides)

    vol_melt_extract_i_func, vol_frac_system_shrinking_tot_func, vol_frac_system_shrinking_func = extract_melt_MAGEMin!(bulk, 0, magemin_interact, out2, grid_model, 1)

    @test bulk == bulk_test
    @test vol_melt_extract_i_func ≈ vol_melt_extract_i
    @test vol_frac_system_shrinking_tot_func ≈ vol_frac_system_shrinking_tot
    @test vol_frac_system_shrinking_func ≈ vol_frac_system_shrinking
end

@testset "MAGEMin crystal fractionation of melt" begin

    data = Initialize_MAGEMin("mp", verbose=false);
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "MnO"; "H2O"; "O"]
    bulk_in_forshaw = [70.48; 25.43/2; 0.77; 3.95; 6.3; 5.54/2; 2.94/2; 0.75/2; 0.07; 100/2; 0.1]

    P = 8.0
    T = 700.0

    sys_in = "mol"

    out     = single_point_minimization(P, T, data, X=bulk_in_forshaw, Xoxides=Xoxides)

    grid_model = CreateGrid(nz=200, Lz=28u"km", tfinal=10.0u"Myr")

    magemin_interact = MAGEMinInteraction(grid=grid_model, crystal_fractionation=true, frac_crystallisation=[0.66, 0.0, 0.0], flux_melting=false, h2o_keep= 0.0, ms_melt_only=false)


    @unpack frac_M, bulk, bulk_M, frac_S, bulk_S, ph, ph_frac, ph_type, ph_id, SS_vec, PP_vec, bulk_wt =  out;
    @unpack depthc_ini, zc = grid_model

    vol_melt_extract_i_func, vol_frac_system_shrinking_tot_func, vol_frac_system_shrinking_func = extract_melt_MAGEMin!(bulk, 0, magemin_interact, out, grid_model, 1)

    Plith = zeros(size(depthc_ini))
    for i in axes(Plith,1)
        Plith[i] = 2700 * 9.81 * reverse(grid_model.zc)[i]
    end

    T0 = zeros(size(depthc_ini))
    for i in axes(zc,1)
        T0[i] = 0.025*zc[i] + 20 + 273.15
    end
    reverse!(T0, dims=1)

    for i in axes(zc,1)
        if reverse(zc)[i] >= 25_000
            T0[i,:] .= 1200 + 273.15
        end
    end


    struct Integrator_Test
        u::Array{Float64, 1}
    end

    integrator = Integrator_Test(T0)

    MeltingCrust.find_ox_index!(magemin_interact.el_index_MAGEMin, magemin_interact.el_index_meltingcrust, magemin_interact.oxides_name, out.oxides)

    @unpack depth_crystal_fractionation, frac_crystallisation = magemin_interact

    bulk_M_ind = @view bulk_M[magemin_interact.el_index_meltingcrust]

    # find index of depth_crystal_fractionation in depth
    depth_index = findmin(abs.(depthc_ini .- depth_crystal_fractionation))[2]

    Plith_frac = round(Plith[depth_index] * 1e-8, digits=2)
    T_frac = round(integrator.u[depth_index] - 273.15, digits=2)

    bulk_M_save = copy(bulk_M)

    # bulk_M_frac = zeros(size(bulk_M))
    # bulk_M_frac[magemin_interact.el_index_MAGEMin] .= bulk_M
    bulk_M_ind .= max.(bulk_M_ind, 0.0001)
    bulk_M_ind .= round.(bulk_M_ind, digits=4)
    # @show bulk_M, Plith_frac, T_frac


    output_melt = single_point_minimization(Plith_frac, T_frac, data, X=bulk_M_ind, Xoxides=Xoxides, sys_in="mol")

    # find all position of "fsp" in output_melt.ph
    fsp_index = findall(x -> x == "fsp", output_melt.ph)
    id_qz = findfirst(output_melt.ph .== "q")
    id_PP_qz = output_melt.ph_id[id_qz] + 1  # +1 because index starts at 0 (C)

    ph_frac_tot = 0  # fraction extracted from the mineral assemblage

    # if there is fsp in the phase, extract composition from output_melt.SS_vec
    if !isempty(fsp_index)
        for j in fsp_index
            # get index for K2O, Na2O and CaO in output_melt.oxides
            k2o_index = findfirst(output_melt.oxides .== "K2O", )
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

    output_melt.bulk .= round.(output_melt.bulk, digits=4)

    bulk_M_new = copy(output_melt.bulk)


    @test sum(bulk_M_save .- bulk_M_new) ≈ 9.71445146547012e-17

    out2     = single_point_minimization(P, T, data, X=bulk_in_forshaw, Xoxides=Xoxides)

    @unpack frac_M, bulk, bulk_M, frac_S, bulk_S, ph, ph_frac, ph_type, ph_id, SS_vec, PP_vec, bulk_wt =  out2;
    @unpack depthc_ini, zc = grid_model

    bulk_M_fonc = copy(bulk_M)
    ph_frac_tot_func = MeltingCrust.process_crystal_fractionation!(bulk_M_fonc, depthc_ini, Plith, integrator, magemin_interact, data, Xoxides)

    @test bulk_M_new == bulk_M_fonc
    @test ph_frac_tot == ph_frac_tot_func

end

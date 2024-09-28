using Unitful
using Parameters

struct CreateGrid
    nz::Int
    Lz::Float64
    Lx::Float64
    Δz::Float64
    Δz_::Float64
    z::Array{Float64, 1}
    zc::Array{Float64, 1}
    depth::Array{Float64, 1}
    depthc_ini::Array{Float64, 1}
    depthc::Array{Float64, 1}
    tfinal::Float64

    function CreateGrid(nz, Lz, Lx, tfinal)
        if Lz <= 0 || Lx <= 0
            error("Length should be positive.")
        elseif tfinal <= 0
            error("Total time should be positive.")
        else
            Δz = Lz / (nz)
            Δz_ = 1 / Δz
            z = collect(range(0, length=nz+1, stop= Lz))
            zc = collect(range(Δz/2, length=nz, stop= Lz-Δz/2))
            depth = reverse(z)
            depthc_ini = reverse(zc)
            depthc = reverse(zc)
            # x first, z second
            new(nz, Lz, Lx, Δz, Δz_, z, zc, depth, depthc_ini, depthc, tfinal)
        end
    end
end

function CreateGrid(;nz::Int, Lz::Unitful.Length, Lx::Unitful.Length= 1u"m", tfinal::Unitful.Time)
    CreateGrid(nz, ustrip(u"m", Lz), ustrip(u"m", Lx), ustrip(u"s", tfinal))
end

@with_kw mutable struct ThermalVariables
    grid::CreateGrid
    convection::BitArray{1} = [true]
    advection::BitArray{1} = [false]
    input_gabbro::BitArray{1} = [true]  # Rayleigh nb
    heat_flux_moho::Float64 = -150*1e-3  # mW.m^-2, very hot Moho
    crust_array::Array{Float64, 1} = zeros(grid.nz)
    neumann_bottom::Bool = true
    time_cooling::Float64 = 2.5e6*365.25*24*3600  # time at which changing bottom boundary
    T0::Array{Float64, 1} = zeros(grid.nz)  # K
    α::Array{Float64, 1} = zeros(grid.nz)  # m^2.s^-1
    Cp::Array{Float64, 1} = zeros(grid.nz)  # J.kg^-1.K^-1
    k::Array{Float64, 1} = zeros(grid.nz)  # W.m^-1.K^-1
    q_conv::Array{Float64, 1} = zeros(grid.nz)  # W.m^-2
    k_conv::Array{Float64, 1} = zeros(grid.nz)  # W.m^-1.K^-1
    Ql::Float64 = 400*1e3  # latent heat kJ/kg
    A::Float64 = 2*1e-6  # radiogenic heat production W.m^-3 (cf Jaupart, C., Mareschal, J.C., 2012. Constraints on crustal heat production from heat flow data, In: Rudnick, R.L. (Ed.) pp. 65–84)
    ρ::Array{Float64, 1} = ones(grid.nz) .* 2_700  # density in kg.m^-3
    Ra::Float64 = (ρ[1] * 9.81 * 3 * 1e-5 * 100 * 8_000^3) / (1e15 * 1e-6)  # Rayleigh nb
    Nu::Array{Float64, 1} = zeros(grid.nz)  # Nusselt nb
end


@with_kw struct PhysicalProperties
    grid::CreateGrid
    Plith::Array{Float64, 1} = zeros(grid.nz)  # lithostatic pressure in Pa
    ϕ::Array{Float64, 1} = zeros(grid.nz)  # porosity
    dϕdT::Array{Float64, 1} = zeros(grid.nz)  # dϕ/dT
    melt_threshold::Float64 = 7  # threshold for melt extraction (in %)
end


@with_kw mutable struct Compositions
    grid::CreateGrid
    # composition of the solid
    compo_solid::Vector{<:Vector{Float64}} = [zeros(11) for _ = 1:grid.nz]
    compo_solid_name::Array{String,1} = fill("", grid.nz)
    compo_solid_matrix::Array{Float64, 2} = zeros(grid.nz, 11)
    # compositions in mol% -> Si, Al, Ca, Mg, Fe, K, Na, Ti, Mn, H, O (thx Jacob)
    compo_pelite_riel::Array{Float64,1} = [63.7,31.43,0.57,4.7,7.08,10.66,4.14,0.67,0.17,100,0]
    # ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "MnO"; "H2O"; "O"]
    compo_pelite_riel_ox::Array{Float64,1} = [compo_pelite_riel[1]; compo_pelite_riel[2]/2; compo_pelite_riel[3]; compo_pelite_riel[4]; compo_pelite_riel[5]; compo_pelite_riel[6]/2; compo_pelite_riel[7]/2; compo_pelite_riel[8]; compo_pelite_riel[9]; compo_pelite_riel[10]/2; compo_pelite_riel[11]]
    compo_semipelite_riel::Array{Float64,1} = [75.21,21.28,0.32,3.58,5.74,5.79,1.7,0.65,0.12,100,0]
    compo_semipelite_riel_ox::Array{Float64,1} = [compo_semipelite_riel[1]; compo_semipelite_riel[2]/2; compo_semipelite_riel[3]; compo_semipelite_riel[4]; compo_semipelite_riel[5]; compo_semipelite_riel[6]/2; compo_semipelite_riel[7]/2; compo_semipelite_riel[8]; compo_semipelite_riel[9]; compo_semipelite_riel[10]/2; compo_semipelite_riel[11]]
    compo_psammite_riel::Array{Float64,1} = [87.33,12.54,0.25,1.33,2.11,2.9,1.61,0.41,0.04,100,0]
    compo_psammite_riel_ox::Array{Float64,1} = [compo_psammite_riel[1]; compo_psammite_riel[2]/2; compo_psammite_riel[3]; compo_psammite_riel[4]; compo_psammite_riel[5]; compo_psammite_riel[6]/2; compo_psammite_riel[7]/2; compo_psammite_riel[8]; compo_psammite_riel[9]; compo_psammite_riel[10]/2; compo_psammite_riel[11]]
    compo_pelite_forshaw::Array{Float64,1} = [70.48,25.43,0.77,3.95,6.3,5.54,2.94,0.75,0.07,100,0]
    compo_pelite_forshaw_ox::Array{Float64,1} = [compo_pelite_forshaw[1]; compo_pelite_forshaw[2]/2; compo_pelite_forshaw[3]; compo_pelite_forshaw[4]; compo_pelite_forshaw[5]; compo_pelite_forshaw[6]/2; compo_pelite_forshaw[7]/2; compo_pelite_forshaw[8]; compo_pelite_forshaw[9]; compo_pelite_forshaw[10]/2; compo_pelite_forshaw[11]]
    # median schists from El Tigre and La Victoria in mol ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "MnO"; "H2O"; "O"]
    compo_pelite_dominguez_ox::Array{Float64,1} = [75.15, 10.58, 0.36, 3.35, 6.04, 3.03, 0.85, 0.6, 0.04, 50.0, 0.302]  # put 10% of Fe2O3 in FeO -> 6.04 FeO becomes 3.02 Fe2O3 so 10% of Fe2O3 is 0.302 O
    compo_semipelite_dominguez_EOET_8_ox::Array{Float64,1} = [81.42, 7.58, 0.20, 2.06, 5.00, 2.20, 0.97, 0.56, 0.01, 50.0, 0.25]
    compo_arkose_pettijohn::Array{Float64,1} = [82.23,14.6,1.51,1.17,2.15,4.65,5.17,0.05,0.05,100,0]
    compo_arkose_pettijohn_ox::Array{Float64,1} = [compo_arkose_pettijohn[1]; compo_arkose_pettijohn[2]/2; compo_arkose_pettijohn[3]; compo_arkose_pettijohn[4]; compo_arkose_pettijohn[5]; compo_arkose_pettijohn[6]/2; compo_arkose_pettijohn[7]/2; compo_arkose_pettijohn[8]; compo_arkose_pettijohn[9]; compo_arkose_pettijohn[10]/2; compo_arkose_pettijohn[11]]
    compo_quartzite_clarke::Array{Float64,1} = [93.09,1.9,0.01,0.26,2.91,2.5,0.14,0.01,0.01,100,0]
    compo_quartzite_clarke_ox::Array{Float64,1} = [compo_quartzite_clarke[1]; compo_quartzite_clarke[2]/2; compo_quartzite_clarke[3]; compo_quartzite_clarke[4]; compo_quartzite_clarke[5]; compo_quartzite_clarke[6]/2; compo_quartzite_clarke[7]/2; compo_quartzite_clarke[8]; compo_quartzite_clarke[9]; compo_quartzite_clarke[10]/2; compo_quartzite_clarke[5] / 2 * 0.1]
    density_pelite::Float64 = 2_500  # kg.m^-3
    density_arkose::Float64 = 2_600  # kg.m^-3
    density_quartzite::Float64 = 2_700  # kg.m^-3
    mass_pelite::Float64 = density_pelite * grid.Lx * grid.Δz    # kg.m^-3
    mass_arkose::Float64 = density_arkose * grid.Lx * grid.Δz    # kg.m^-3
    mass_quartzite::Float64 = density_quartzite * grid.Lx * grid.Δz  # kg.m^-3
    mass_frac_melt::Array{Float64,1} = zeros(grid.nz)  # mass fraction of melt
end


@with_kw mutable struct MAGEMinInteraction
    grid::CreateGrid
    miniz_calls::Array{Float64, 1} = zeros(grid.nz)
    updated::Array{Bool, 1} = zeros(grid.nz)
    ΔT::Array{Float64, 1} = zeros(grid.nz)
    depth_min_thermo::Float64 = 8_500  # minimum depth for thermo
    depth_max_thermo::Float64 = 25_000  # maximum depth for thermo
    mineral_content::Array{String,1} = [""]
    oxides_name::Array{String,1} = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "MnO"; "H2O"; "O"]
    el_index_MAGEMin::Array{Int, 1} = zeros(size(oxides_name,1))
    el_index_meltingcrust::Array{Int, 1} = zeros(size(oxides_name,1))
    vol_melt_extract::Array{Float64, 2} = zeros(grid.nz, 1)
    vol_frac_system_shrinking_tot::Array{Float64, 2} = zeros(grid.nz, 1)
    compo_melt_wt::Array{Float64, 2} = zeros(grid.nz, size(oxides_name,1))  # composition of the melt
    compo_melt_wt_all::Vector{<:Vector{Float64}} = Vector{Float64}[] # composition of the melt
    compo_melt_wt_all_index::Vector{Int} = Int[] # composition of the melt
    compo_melt_all_time::Vector{Float64} = Float64[] # composition of the melt
    compo_melt_all_volume::Vector{Float64} = Float64[] # volume of the melt
    compo_melt_all_unit::Vector{String} = String[] # unit of the melt
    compo_melt_i_wt::Array{Float64, 1} = zeros(size(oxides_name,1))  # composition of the melt
    compo_melt_wt_dry_i::Array{Float64, 1} = zeros(size(oxides_name,1)-2)  # composition of the melt
    compo_melt_wt_dry_all::Vector{<:Vector{Float64}} = Vector{Float64}[]  # composition of the melt
    compo_i_g::Array{Float64, 1} = zeros(size(oxides_name,1))  # composition of the melt
    compo_rock_wt::Array{Float64, 2} = zeros(grid.nz, size(oxides_name,1))  # composition of the melt
    time_melt::Array{Float64, 1} = zeros(1)  # array to save the time of melt extraction
    time_volume::Array{Float64, 1} = zeros(1)  # array to save the time of melt extraction
    total_volume_melt::Float64 = 0  # total volume of melt extracted
    name_assemblage::Matrix{String} = fill("", grid.nz, 11)  # name of the assemblage
    mode_assemblage::Matrix{Float64} = zeros(grid.nz, 11)  # name of the assemblage
    mode_assemblage_vol::Matrix{Float64} = zeros(grid.nz, 11)  # name of the assemblage
    compo_assemblage::Matrix{Float64} = zeros(length(oxides_name),grid.nz*11)  # name of the assemblage
    crystal_fractionation::Bool = false  # crystal fractionation
    depth_crystal_fractionation::Float64 = 13_000  # depth for crystal fractionation (m)
    frac_crystallisation::Array{Float64, 1} = [0, 0, 0]  # fraction of crystallisation (pl, kfs, qz)
    flux_melting::Bool = false  # if flux melting is activated or not
    melt_prev::Array{Bool, 1} = zeros(grid.nz) # if melt was extracted at the previous timestep or not
    h2o_keep::Float64 = 0  # keep 0% of water in the bulk rock for water flux melting
    ms_melt_only::Bool = false
    ms_melt_check::Array{Bool, 1} = ones(grid.nz)  # if melt and ms are present
    grt_frac::Array{Float64, 1} = zeros(grid.nz)  # fraction of garnet
    grt_frac_prev::Array{Float64, 1} = zeros(grid.nz)  # fraction of garnet at the previous timestep
    grt_frac_remove::Float64 = 0.0  # fraction of garnet to remove
    crd_frac::Array{Float64, 1} = zeros(grid.nz)  # fraction of cordierite
    crd_frac_prev::Array{Float64, 1} = zeros(grid.nz)  # fraction of cordierite at the previous timestep
    crd_frac_remove::Float64 = 0.0  # fraction of cordierite to remove
end
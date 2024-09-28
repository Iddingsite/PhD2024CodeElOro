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
    M::Float64 = 1  # melt-focusing factor
    crust_array::Array{Float64, 1} = zeros(grid.nz)
    neumann_bottom::Bool = true
    melt_emplaced::Float64 = 0  # melt emplaced in the crust
    time_cooling::Float64 = 2.5e6*365.25*24*3600  # time at which changing bottom boundary
    fertile_lithosphere::Array{Bool, 1} = fill(true, grid.nz)
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


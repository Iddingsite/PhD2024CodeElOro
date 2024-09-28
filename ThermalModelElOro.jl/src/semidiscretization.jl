
function stencil_diffusion!(dtT, T, k, Cp, ρ, A, z, input_gabbro, heat_flux_moho, neumann_bottom)

    @inline qz(D,T,iz,z) = (0.5*(D[iz]+D[iz+1])) * (T[iz+1]-T[iz]) / (z[iz+1]-z[iz])

    # iterate inside the arrays
    # Threads.@threads for iz=2:size(dtT,1)-1
    for iz=2:size(dtT,1)-1
        dtT[iz] = 1/(ρ[iz]*Cp[iz]) *
                  ((qz(k,T,iz,z) - qz(k,T,iz-1,z)) / (z[iz+1]-z[iz]) + A)
    end

    dtT[end] = 0


    if neumann_bottom == true
        if input_gabbro[1] == false
            # neumann bottom boundary
            dtT[1] = 1/(ρ[1]*Cp[1]) * (((qz(k,T,1,z) - heat_flux_moho)/ (z[2]-z[1])) + A)
        end
    else
        dtT[1] = 0
    end

    return
end

function semi_discretisation(du,u,p,t)

    @unpack geoparam, grid, thermal_param, physical_param = p
    @unpack z, depthc = grid
    @unpack α, k, Cp, ρ, Ql, A, k_conv, q_conv, Nu, convection, advection, input_gabbro, heat_flux_moho, neumann_bottom = thermal_param
    @unpack ϕ, dϕdT = physical_param
    MatParam, phases = geoparam

    # update Cp
    compute_heatcapacity!(Cp, MatParam, phases, (;T=u))
    compute_conductivity!(k, MatParam, phases, (;T=u))

    # p_rhyolite = MeltingParam_Smooth3rdOrder(a=3043.0, b=-10552.0, c=12204.9, d= -4709.0)
    p_rhyolite = MeltingParam_Quadratic(T_s=(680 + 273.15), T_l=(1150 + 273.15))

    # p_basalt = MeltingParam_Smooth3rdOrder(a=517.9, b=-1619.0, c=1699.0, d=-597.4)
    p_basalt = MeltingParam_Quadratic(T_s=(820 + 273.15), T_l=(1200 + 273.15))

    # find last index of bellow 30_000 in depthc
    boundary = findall(depthc .>= 30_000)


    # create views

    # ϕ_basalt = @view ϕ[2:boundary]
    # dϕdT_basalt = @view dϕdT[2:boundary]
    # ϕ_rhyolite = @view ϕ[boundary+1:end]
    # dϕdT_rhyolite = @view dϕdT[boundary+1:end]
    # T_bottom = @view u[2:boundary]
    # T_top = @view u[boundary+1:end]

    # # calculate the melt fraction
    # compute_meltfraction!(ϕ_basalt, p_basalt, (;T=T_bottom))
    # compute_dϕdT!(dϕdT_basalt, p_basalt, (;T=T_bottom))
    # compute_meltfraction!(ϕ_rhyolite, p_rhyolite, (;T=T_top))
    # compute_dϕdT!(dϕdT_rhyolite, p_rhyolite, (;T=T_top))

    # check if boundary is empty or not
    # if it is, use only rhyolite, if it is not use both

    if isempty(boundary)
        # calculate the melt fraction
        compute_meltfraction!(ϕ, p_rhyolite, (;T=u))
        compute_dϕdT!(dϕdT, p_rhyolite, (;T=u))
    else
        ϕ_basalt = @view ϕ[2:boundary[1]]
        dϕdT_basalt = @view dϕdT[2:boundary[1]]
        ϕ_rhyolite = @view ϕ[boundary[1]+1:end]
        dϕdT_rhyolite = @view dϕdT[boundary[1]+1:end]
        T_bottom = @view u[2:boundary[1]]
        T_top = @view u[boundary[1]+1:end]

        # calculate the melt fraction
        compute_meltfraction!(ϕ_basalt, p_basalt, (;T=T_bottom))
        compute_dϕdT!(dϕdT_basalt, p_basalt, (;T=T_bottom))
        compute_meltfraction!(ϕ_rhyolite, p_rhyolite, (;T=T_top))
        compute_dϕdT!(dϕdT_rhyolite, p_rhyolite, (;T=T_top))
    end

    # implement latent heat
    Cp .= Cp .+ (dϕdT .* Ql)

    # semi-discretization
    if convection[1] == true
        convective_heat_flux!(q_conv, Nu, k)
        convective_thermal_conductivity!(k_conv, q_conv, k)  # W.m^-1.K^-1
        stencil_diffusion!(du, u, k_conv, Cp, ρ, A, z, input_gabbro, heat_flux_moho, neumann_bottom)
    elseif convection[1] == false
        if advection[1] == true
            stencil_diffusion!(du, u, k, Cp, ρ, A, z, input_gabbro, heat_flux_moho, neumann_bottom)
        elseif advection[1] == false
            stencil_diffusion!(du, u, k, Cp, ρ, A, z, input_gabbro, heat_flux_moho, neumann_bottom)
        end
    else
        error("wrong settings")
    end

end


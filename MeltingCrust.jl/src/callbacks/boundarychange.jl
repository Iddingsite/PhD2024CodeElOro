using Parameters

function change_boundary!(integrator)

    @unpack thermal_param = integrator.p
    @unpack input_gabbro = thermal_param

    @show input_gabbro
    input_gabbro[1] = 0
    @show input_gabbro

end

function change_boundary_convection_gabbro!(integrator)

    @unpack thermal_param, grid = integrator.p
    @unpack input_gabbro = thermal_param

    # remove neumann boundary condition
    thermal_param.neumann_bottom = false
    # emplace gabbro
    boundary_top = findall(grid.depthc .>= 25_000)[end]
    integrator.u[1:boundary_top] .= 1250 + 273.15
    # activate convection
    thermal_param.convection[1] = true

end

function change_boundary_convection_blueschist!(integrator)

    @unpack thermal_param = integrator.p
    @unpack input_gabbro = thermal_param

    # emplace blueschist
    integrator.u[1] = 300 + 273.15
    # activate convection
    thermal_param.convection[1] = false
    thermal_param.neumann_bottom = false

    println("hi")

    sleep(10)
end

function change_boundary_convection_flux!(u,t,integrator)

    @unpack thermal_param = integrator.p
    @unpack heat_flux_moho = thermal_param

    # convert time from s to Myr
    t_Ma = t/(3600*24*365.25*1e6)

    if t_Ma <= 8
        thermal_param.heat_flux_moho = -(7*t_Ma + 80)*1e-3
    end

end
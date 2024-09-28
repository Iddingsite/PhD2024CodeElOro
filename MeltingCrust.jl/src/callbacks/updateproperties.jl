
# update thermal properties for explicit schemes
function thermal_param_call_func(u, t, integrator)

    @unpack grid, thermal_param = integrator.p
    @unpack α, Cp, k, q_conv, k_conv, ρ, Nu = thermal_param

    thermal_diffusivity!(α, u)  # m^2.s^-1
    heat_capacity!(Cp, u)  # J.kg^-1.K^-1
    thermal_conductivity!(k, α, Cp, ρ)  # W.m^-1.K^-1
    convective_heat_flux!(q_conv, Nu, k)
    convective_thermal_conductivity!(k_conv, q_conv, k)  # W.m^-1.K^-1
end

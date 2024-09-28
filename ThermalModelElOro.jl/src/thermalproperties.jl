"
Thermal conductivity
Calculated from Whittington et al, 2009 in m2/s
"
function thermal_diffusivity!(α, T)
    for iz=eachindex(α)
        if T[iz] < 846
            α[iz] = (567.3 / T[iz] - 0.062) * 1e-6
        else
            α[iz] = (0.732 - 0.000135 * T[iz]) * 1e-6
        end
    end
    return
end


"
Specific heat capacity
Calculated from Whittington et al, 2009 in J.kg^-1.K^-1
"
function heat_capacity!(Cp, T)
    for iz=eachindex(Cp)
        if T[iz] < 846
            Cp[iz] = (199.50 + 0.0857 * T[iz] - 5.0*1e6 * T[iz]^(-2)) ./ 221.78 * 1000
        else
            Cp[iz] = (229.32 + 0.0323 * T[iz] - 47.9*1e-6 * T[iz]^(-2)) ./ 221.78 * 1000
        end
    end
    return
end

"
Specific heat capacity for lithosphere
Calculated from McKenzie and al., 2005, assuming 0.11 mol% of fayalite and 0.89 mol% of forsterite
"
function heat_capacity_lithosphere!(Cp, T)
    for iz=eachindex(Cp)
        Cp[iz] = 233
    end
    return
end



"
Thermal conductivity
Calculated from McKenzie and al., 2005 in m2/s
"
function thermal_conductivity_lithosphere!(k, T)
    for iz=eachindex(k)
        k[iz] = 5.3 / (1 + 0.0015 * (T[iz]-273.15)) + (1.753*1e-2) / (0+1) * (T[iz]-273.15 + 273.15)^(0+1) +
                -1.0365*1e-4 / (1+1) * (T[iz]-273.15 + 273.15)^(1+1) + 2.2451*1e-7 / (2+1) * (T[iz]-273.15 + 273.15)^(2+1) + -3.4071*1e-11 / (3+1) * (T[iz]-273.15 + 273.15)^(3+1)
    end
    return
end


"
Thermal conductivity
Calculated from Whittington et al, 2009 in m2/s
"
function thermal_conductivity_crust!(k, α, Cp, ρ)
    k = α * ρ * Cp
    return k
end


"
Thermal conductivity
Calculated from Whittington et al, 2009 in m2/s
"
function thermal_conductivity!(k, α, Cp, ρ)
    for iz=eachindex(k)
        k[iz] = α[iz] * ρ[iz] * Cp[iz]
    end
    return
end


"
Convective heat flux
"
function convective_heat_flux!(q_conv, Nu, k)
    for iz=eachindex(q_conv)
        q_conv[iz] = Nu[iz] * k[iz] * 100 / 8_000
    end
    return
end


"
thermal conductivity including convective heat flux
kconv = h*Q/ΔT + k
"
function convective_thermal_conductivity!(k_conv, q_conv, k)
    for iz=eachindex(k_conv)
        k_conv[iz] = 8000 * q_conv[iz] / 100 + k[iz]
    end
    return
end
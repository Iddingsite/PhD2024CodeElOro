using Random
# define seed
Random.seed!(1234)
using Distributions


function remove_and_emplace_melt_func(u,t,integrator)

    @unpack grid, physical_param, thermal_param = integrator.p
    @unpack depth, depthc, depthc_ini, nz, z, zc = grid
    @unpack Plith, ϕ = physical_param
    @unpack M, melt_emplaced = thermal_param

    # if t < ustrip(u"s",20u"Myr")
    # p_rhyolite = MeltingParam_Smooth3rdOrder(a=3043.0, b=-10552.0, c=12204.9, d= -4709.0)
    p_rhyolite = MeltingParam_Quadratic(T_s=(680 + 273.15), T_l=(1150 + 273.15))
    # p_basalt = MeltingParam_Smooth3rdOrder(a=517.9, b=-1619.0, c=1699.0, d=-597.4)
    p_basalt = MeltingParam_Quadratic(T_s=(820 + 273.15), T_l=(1200 + 273.15))

    # find last index of bellow 30_000 in depthc
    # boundary = findall(depthc .>= 30_000)[end]
    boundary = 1

    # create views
    ϕ_basalt = @view ϕ[1:boundary]
    ϕ_rhyolite = @view ϕ[boundary+1:end]
    T_bottom = @view u[1:boundary]
    T_top = @view u[boundary+1:end]

    # calculate the melt fraction
    compute_meltfraction!(ϕ_basalt, p_basalt, (;T=T_bottom))
    compute_meltfraction!(ϕ_rhyolite, p_rhyolite, (;T=T_top))

    # compute_meltfraction!(ϕ, p_rhyolite, (;T=u))

    # if ϕ is higher than 7%, extract the melt and add it at 10 km above the point
    # mix the temperature of the melt with the solid at the depth
    # shift the column of rocks bellow the point

    # loop over the melt fraction
    # for i in axes(ϕ,1)
    i = 1
    if ϕ[i] > 0.07
        # extract the melt
        melt = ϕ[i]
        ϕ[i] = 0

        # find the closest depth to depth_current - 10000
        # index_emplacement = argmin(abs.(depthc .- (depth_current - 10000)))
        # index_emplacement = argmin(abs.(depthc .- (20000)))

        # # define randomly a depth between 180000 and 20000
        # index_emplacement = rand(17000:19000)

        # define a normal distribution around 18000 with a standard deviation of 1000
        depth_emplacement = rand(Normal(19000, 1500))
        
        # depth_emplacement = 22000

        index_emplacement = argmin(abs.(depthc .- (depth_emplacement)))

        # mix the temperature of the melt with the solid at the depth
        T_melt = u[i]
        T_solid = u[index_emplacement]

        # mix the temperature with relative volume of melt and rock
        dz_melt = z[i+1]-z[i]
        dz_hostrock = z[index_emplacement+1]-z[index_emplacement]

        u[index_emplacement] = (T_melt * melt * dz_melt * M + T_solid * 1 * dz_hostrock) / (melt * dz_melt * M + 1 * dz_hostrock)

        thermal_param.melt_emplaced += (melt * dz_melt)

        # shift the Temperature of the column of rocks bellow the point by melt * z[i+1]-z[i] length using advection equation
        # u_new = u_old - v * dt * (u_old - u_old-1) / dz
        # where v is the velocity of the column. Here we assume the velocity is ϕ[i] * dz / dt where dt is 1 (basically intantenious transport)
        # and u_old-1 is u[i-1]

        u_old = copy(u)

        for j in 2:index_emplacement
            u[j] = u_old[j] - melt * dz_melt * (u_old[j] - u_old[j-1]) / (z[j] - z[j-1])
        end


        # replace temperature of the point of extraction to solidus temperature
        # check if it is bellow or above boundary
        if i <= boundary
            u[i] = 820 + 273.15
        else
            u[i] = 770 + 273.15
        end

    end
    # end


    # end

    # display(plot(zc, u .- 273.15, label="Temperature (°C)", xlabel="Depth (m)", ylabel="Temperature (°C)", title="Temperature profile", legend=:topleft))

    # # sleep for 0.1 seconds
    # sleep(0.1)
end

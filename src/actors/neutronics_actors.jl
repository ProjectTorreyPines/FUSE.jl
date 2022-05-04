mutable struct neutron_particle
    x
    y
    z
    δvx
    δvy
    δvz
end

function Rcoord(n::neutron_particle)
    sqrt(n.x^2 + n.y^2)
end

function Zcoord(n::neutron_particle)
    n.z
end

#= =============== =#
#  ActorNeutronics  #
#= =============== =#

mutable struct ActorNeutronics <: ActorAbstract
    dd::IMAS.dd
end

function ParametersActor(::Type{Val{:ActorNeutronics}})
    par = ParametersActor(nothing)
    par.N = Entry(Integer, "", "Number of particles"; default=100000)
    par.step = Entry(Real, "", "Interator stepping"; default=0.05)
    par.do_plot = Entry(Bool, "", "plot"; default=false)
    return par
end

function ActorNeutronics(dd::IMAS.dd, act::ParametersActor; kw...)
    par = act.ActorNeutronics(kw...)
    actor = ActorNeutronics(dd)
    step(actor; N=par.N, step=par.step, do_plot=par.do_plot)
    finalize(actor)
    return actor
end

function step(actor::ActorNeutronics; N::Integer=100000, step=0.05, do_plot::Bool=false)
    dd = actor.dd
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]

    # 2D neutron source
    neutron_source_1d = IMAS.alpha_heating(cp1d) * 4
    neutron_source_2d = transpose(IMAS.interp1d(cp1d.grid.psi, neutron_source_1d).(eqt.profiles_2d[1].psi))

    # 2D neutron source (masked)
    r = eqt.profiles_2d[1].grid.dim1
    z = eqt.profiles_2d[1].grid.dim2
    rz_lcfs = collect(zip(eqt.boundary.outline.r, eqt.boundary.outline.z))
    mask = [IMAS.PolygonOps.inpolygon((rr, zz), rz_lcfs) for zz in z, rr in r]
    R = [rr for zz in z, rr in r]
    Z = [zz for zz in z, rr in r]
    neutron_source_2d .= neutron_source_2d .* mask

    # cumulative distribution function
    CDF = cumsum(neutron_source_2d[:])
    count_per_particle = CDF[end] / N
    CDF .= (CDF .- CDF[1]) ./ (CDF[end] - CDF[1])
    draw = Int.(ceil.(IMAS.interp1d(CDF, 1:length(CDF), :linear).(rand(N))))

    # neutron structures
    ϕ = rand(N) .* 2pi
    θv = rand(N) .* 2pi
    ϕv = acos.(rand(N) .* 2.0 .- 1.0)
    neutrons = neutron_particle.(
        R[draw] .* cos.(ϕ),
        R[draw] .* sin.(ϕ),
        Z[draw],
        step * sin.(ϕv) .* cos.(θv),
        step * sin.(ϕv) .* sin.(θv),
        step * cos.(ϕv)
    )

    wall = IMAS.first_wall(dd.wall)

    # plot birth of neutrons
    if do_plot
        histogram2d(
            Rcoord.(neutrons),
            Zcoord.(neutrons),
            nbins=(LinRange(minimum(r), maximum(r), length(r) - 1), LinRange(minimum(z), maximum(z), length(z) - 1)),
            aspect_ratio=:equal)
        display(plot!(wall.r, wall.z, label="", title="Neutron start"))
    end

    # advance neutrons until they hit the wall
    rz_wall = collect(zip(wall.r, wall.z))
    Threads.@threads for n in neutrons
        while FUSE.PolygonOps.inpolygon((Rcoord(n), Zcoord(n)), rz_wall) == 1
            n.x += n.δvx
            n.y += n.δvy
            n.z += n.δvz
        end
    end

    # plot neutrons hit the wall
    if do_plot
        histogram2d(
            Rcoord.(neutrons),
            Zcoord.(neutrons),
            nbins=(LinRange(minimum(r), maximum(r), length(r) - 1), LinRange(minimum(z), maximum(z), length(z) - 1)),
            aspect_ratio=:equal)
        display(plot!(wall.r, wall.z, label="", title="Neutron stop"))
    end

    # find counts/s/m²
    # smooth the load of each neutron withing a window
    # `window` approach only works only for equispaced wall segments
    wall_r, wall_z = IMAS.resample_2d_line(wall.r, wall.z)
    wall_r, wall_z = wall_r[1:end-1], wall_z[1:end-1]
    d = sqrt.(IMAS.diff(vcat(wall_r, wall_r[1])) .^ 2.0 .+ IMAS.diff(vcat(wall_z, wall_z[1])) .^ 2.0)
    d = (d + vcat(d[end], d[1:end-1])) / 2.0
    s = 1.0 ./ (d .* wall_r .* 2pi)
    stencil = collect(-10:10)
    window = 1.0 ./ (abs.(stencil) .+ 1.0) .- (1.0 / (maximum(stencil) + 1))
    window = window ./ sum(window) .* count_per_particle
    N_sm² = zeros(size(wall_r))
    for n in neutrons
        index = argmin((wall_r .- Rcoord(n)) .^ 2.0 + (wall_z .- n.z) .^ 2.0)
        index = mod.(stencil .+ index .- 1, length(wall_r)) .+ 1
        tmp = s[index] ./ sum(s[index]) .* window
        N_sm²[index] .+= tmp
    end

    # IMAS assignements
    dd.neutronics.first_wall.r = wall_r
    dd.neutronics.first_wall.z = wall_z
    resize!(dd.neutronics.time_slice, 1)
    dd.neutronics.time_slice[1].wall_loading = N_sm²

    # # plot neutron wall loading
    # if do_plot
    #     display(plot(wall_r, wall_z, line_z=N_sm², aspect_ratio=:equal, linewidth=10, label="", xlim=[minimum(wall_r) * 0.9, maximum(wall_r) * 1.1], clim=(0.0, maximum(N_sm²))))
    # end

    # # plot neutron wall loading
    # if do_plot
    #     display(plot(cumsum(d) / sum(d), N_sm², label="Neutron wall loading"))
    # end

    # # wall displacement based on neutrons' incoming direction
    # wall_dr0 = zeros(size(wall_r))
    # wall_dz0 = zeros(size(wall_z))
    # for n in neutrons
    #     index = argmin((wall_r .- Rcoord(n)) .^ 2.0 + (wall_z .- Zcoord(n)) .^ 2.0)
    #     old_r = Rcoord(n)
    #     old_z = Zcoord(n)
    #     n.x += n.δvx
    #     n.y += n.δvy
    #     n.z += n.δvz
    #     wall_dr0[index] += (Rcoord(n) - old_r)
    #     wall_dz0[index] += (Zcoord(n) - old_z)
    # end
    # wall_r0 = wall_r .+ 0.0
    # wall_z0 = wall_z .+ 0.0

    # # smooth
    # wall_dr0 = IMAS.smooth(wall_dr0 .* scale, length(wall_dr0) / 4)
    # wall_dz0 = IMAS.smooth(wall_dz0 .* scale, length(wall_dr0) / 4)

    # # scale
    # scale = diff(collect(extrema(wall_r0))) ./ diff(collect(extrema(wall_dr0))) / 1
    # wall_dr0 .*= scale
    # wall_dz0 .*= scale

    # # displace
    # wall_r = vcat(wall_r0, wall_r0[1])
    # wall_dr = vcat(wall_dr0, wall_dr0[1])
    # wall_z = vcat(wall_z0, wall_z0[1])
    # wall_dz = vcat(wall_dz0, wall_dz0[1])

    # if do_plot
    #     plot(wall_r, wall_z, aspect_ratio=:equal, label="First wall")
    #     plot!(wall_r .+ wall_dr, wall_z .+ wall_dz, label="Blanket")
    #     plot!(xlim=[minimum(wall_r) * 0.5, maximum(wall_r) * 1.5])
    # end



end
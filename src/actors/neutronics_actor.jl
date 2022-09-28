#= =============== =#
#  ActorNeutronics  #
#= =============== =#

mutable struct neutron_particle{T<:Real}
    x::T
    y::T
    z::T
    δvx::T
    δvy::T
    δvz::T
end

function Rcoord(n::neutron_particle)
    sqrt(n.x^2 + n.y^2)
end

function Zcoord(n::neutron_particle)
    n.z
end

mutable struct ActorNeutronics <: PlasmaAbstractActor
    dd::IMAS.dd
    par::ParametersActor
    function ActorNeutronics(dd::IMAS.dd, par::ParametersActor; kw...)
        logging_actor_init(ActorNeutronics)
        par = par(kw...)
        return new(dd, par)
    end
end

function ParametersActor(::Type{Val{:ActorNeutronics}})
    par = ParametersActor(nothing)
    par.N = Entry(Integer, "", "Number of particles"; default=100000)
    par.step = Entry(Float64, "", "Interator stepping"; default=0.05)
    par.do_plot = Entry(Bool, "", "plot"; default=false)
    return par
end

"""
    ActorNeutronics(dd::IMAS.dd, act::ParametersAllActors; kw...)

This actor estimates the neutron loading on the wall using the fusion source from `dd.core_sources`.

!!! note 
    Stores data in `dd.neutronics`
"""
function ActorNeutronics(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorNeutronics(kw...)
    actor = ActorNeutronics(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorNeutronics; N::Integer=actor.par.N, step::Float64=actor.par.step, do_plot::Bool=actor.par.do_plot)
    dd = actor.dd
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]

    # 2D neutron source
    neutron_source_1d = IMAS.alpha_heating(cp1d) * 4 # W/m^3
    neutron_source_2d = transpose(IMAS.interp1d(cp1d.grid.psi, neutron_source_1d).(eqt.profiles_2d[1].psi)) # W/m^3

    # 2D neutron source (masked)
    r = eqt.profiles_2d[1].grid.dim1
    z = eqt.profiles_2d[1].grid.dim2
    rz_lcfs = collect(zip(eqt.boundary.outline.r, eqt.boundary.outline.z))
    mask = [IMAS.PolygonOps.inpolygon((rr, zz), rz_lcfs) for zz in z, rr in r]
    R = [rr for zz in z, rr in r]
    Z = [zz for zz in z, rr in r]
    neutron_source_2d .= neutron_source_2d .* mask
    neutron_source_2d .*= R .* (2pi * (r[2] - r[1]) * (z[2] - z[1])) # W

    # cumulative distribution function
    CDF = cumsum(neutron_source_2d[:])
    W_per_trace = CDF[end] / N
    CDF .= (CDF .- CDF[1]) ./ (CDF[end] - CDF[1])
    draw = Int.(ceil.(IMAS.interp1d(CDF, 1:length(CDF), :linear).(rand(N))))

    # neutron structures
    ϕ = rand(N) .* 2pi
    θv = rand(N) .* 2pi
    ϕv = acos.(rand(N) .* 2.0 .- 1.0)
    neutrons = neutron_particle.(R[draw] .* cos.(ϕ), R[draw] .* sin.(ϕ), Z[draw], step * sin.(ϕv) .* cos.(θv), step * sin.(ϕv) .* sin.(θv), step * cos.(ϕv))

    # resample wall and make sure it's clockwise (for COCOS = 11)
    wall = IMAS.first_wall(dd.wall)
    wall_r, wall_z = IMAS.resample_2d_line(wall.r, wall.z; step=0.1)
    R0 = eqt.global_quantities.magnetic_axis.r
    Z0 = eqt.global_quantities.magnetic_axis.z
    IMAS.reorder_flux_surface!(wall_r, wall_z, R0, Z0)

    # advance neutrons until they hit the wall
    rz_wall = collect(zip(wall_r, wall_z))
    Threads.@threads for n in neutrons
        while PolygonOps.inpolygon((Rcoord(n), Zcoord(n)), rz_wall) == 1
            n.x += n.δvx
            n.y += n.δvy
            n.z += n.δvz
        end
    end

    # find neutron flux [counts/s/m²]
    # smooth the load of each neutron withing a window
    wall_r, wall_z = wall_r[1:end-1], wall_z[1:end-1]
    d = sqrt.(IMAS.gradient(vcat(wall_r, wall_r[1])) .^ 2.0 .+ IMAS.gradient(vcat(wall_z, wall_z[1])) .^ 2.0)
    d = (d[1:end-1] .+ d[2:end]) / 2.0
    l = cumsum(d)
    wall_s = d .* wall_r .* 2π
    ns = 10
    stencil = collect(-ns:ns)

    nflux_r = zeros(size(wall_r))
    nflux_z = zeros(size(wall_z))
    ll = stencil .* 0.0
    for n in neutrons
        old_r = Rcoord(n)
        old_z = Zcoord(n)
        index0 = argmin((wall_r .- old_r) .^ 2.0 + (wall_z .- old_z) .^ 2.0)
        index = mod.(stencil .+ index0 .- 1, length(wall_r)) .+ 1

        n.x += n.δvx
        n.y += n.δvy
        n.z += n.δvz
        new_r = Rcoord(n)
        new_z = Zcoord(n)

        ll .= cumsum(d[index])
        ll .-= ll[ns+1]
        window = exp.(-(ll ./ (l[end] / length(l)) / (2 * ns + 1) * 5) .^ 2)
        window = window ./ sum(window)
        unit_vector = sqrt((new_r - old_r)^2 + (new_z - old_z)^2)

        nflux_r[index] += (new_r - old_r) ./ unit_vector .* window .* W_per_trace ./ wall_s[index]
        nflux_z[index] += (new_z - old_z) ./ unit_vector .* window .* W_per_trace ./ wall_s[index]
    end

    # IMAS assignements
    dd.neutronics.first_wall.r = wall_r
    dd.neutronics.first_wall.z = wall_z
    ntt = resize!(dd.neutronics.time_slice)
    ntt.wall_loading.flux_r = nflux_r
    ntt.wall_loading.flux_z = nflux_z
    ntt.wall_loading.power = sqrt.(nflux_r .^ 2.0 .+ nflux_z .^ 2.0) .* wall_s

    # plot neutron wall loading cx
    if do_plot
        neutrons = neutron_particle.(R[draw] .* cos.(ϕ), R[draw] .* sin.(ϕ), Z[draw], step * sin.(ϕv) .* cos.(θv), step * sin.(ϕv) .* sin.(θv), step * cos.(ϕv))

        histogram2d(
            Rcoord.(neutrons),
            Zcoord.(neutrons),
            nbins=(LinRange(minimum(r), maximum(r), length(r) - 1), LinRange(minimum(z), maximum(z), length(z) - 1)),
            aspect_ratio=:equal,
            weights=zeros(N) .+ 1 / 40,
        )

        plot!(ntt.wall_loading, xlim=[minimum(wall.r) * 0.9, maximum(wall.r) * 1.1], title="First wall neutron loading")
        display(plot!(wall.r, wall.z, color=:black, label="", xlim=[minimum(wall.r) * 0.9, maximum(wall.r) * 1.1], title=""))
    end

end
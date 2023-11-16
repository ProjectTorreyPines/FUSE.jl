#= =============== =#
#  ActorNeutronics  #
#= =============== =#
Base.@kwdef mutable struct FUSEparameters__ActorNeutronics{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    N::Entry{Int} = Entry{Int}("-", "Number of particles"; default=100000)
    do_plot::Entry{Bool} = Entry{Bool}("-", "Plot"; default=false)
end

mutable struct ActorNeutronics{D,P} <: PlasmaAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorNeutronics{P}
    function ActorNeutronics(dd::IMAS.dd{D}, par::FUSEparameters__ActorNeutronics{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorNeutronics)
        par = par(kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorNeutronics(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the neutron wall loading

!!! note

    Stores data in `dd.neutronics`
"""
function ActorNeutronics(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorNeutronics(dd, act.ActorNeutronics; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorNeutronics)
    do_plot::Bool = actor.par.do_plot
    p = do_plot ? plot() : nothing

    neutrons, W_per_trace = define_neutrons(actor)
    if do_plot
        plot!(p, neutrons, actor.dd.equilibrium.time_slice[])
    end

    rwall, zwall = define_wall(actor)
    rz_wall = collect(zip(rwall, zwall))
    Nw = length(rz_wall)

    # advance neutrons until they hit the wall
    for n in neutrons
        #println(n)
        ti = Inf

        xn = n.x
        yn = n.y
        zn = n.z
        vx = n.δvx
        vy = n.δvy
        vz = n.δvz

        v2 = vx^2 + vy^2
        vz2 = vz^2

        a = v2
        b = 2 * (vx * xn + vy * yn)
        c = xn^2 + yn^2
        rmin = sqrt(c - b^2 / (4a))

        @inbounds for k in eachindex(rz_wall)
            rw1, zw1 = rz_wall[k]
            rw2, zw2 = k == Nw ? rz_wall[1] : rz_wall[k+1]

            rw1 == rw2 && zw1 == zw2 && continue

            if zw1 > zw2
                rw1, rw2 = rw2, rw1
                zw1, zw2 = zw2, zw1
            end

            (vz > 0 && zn > zw2) && continue
            (vz < 0 && zn < zw1) && continue
            (rw1 < rmin && rw2 < rmin) && continue

            t = toroidal_intersection(rw1, zw1, rw2, zw2, xn, yn, zn, vx, vy, vz, v2, vz2)
            if t < ti
                ti = t
            end
        end
        n.x += vx * ti
        n.y += vy * ti
        n.z += vz * ti
    end

    # find neutron flux [counts/s/m²]
    # smooth the load of each neutron withing a window
    # note: flux is defined at the cells, not at the nodes
    wall_r = (rwall[1:end-1] .+ rwall[2:end]) ./ 2.0
    wall_z = (zwall[1:end-1] .+ zwall[2:end]) ./ 2.0
    d = sqrt.(IMAS.gradient(rwall) .^ 2.0 .+ IMAS.gradient(zwall) .^ 2.0)
    d = @views (d[1:end-1] .+ d[2:end]) ./ 2.0
    l = cumsum(d)
    wall_s = d .* wall_r .* 2π
    ns = 10 # window size
    stencil = collect(-ns:ns)

    nflux_r = zero(wall_r)
    nflux_z = zero(wall_z)
    dwall = zero(wall_r)
    index = zero(stencil)
    ll = zeros(2ns + 1)
    window = zeros(2ns + 1)
    for n in neutrons
        old_r = Rcoord(n)
        old_z = Zcoord(n)
        @. dwall = (wall_r - old_r)^2 + (wall_z - old_z)^2
        index0 = argmin(dwall)
        index .= mod.(stencil .+ index0 .- 1, length(wall_r)) .+ 1

        n.x += n.δvx
        n.y += n.δvy
        n.z += n.δvz
        new_r = Rcoord(n)
        new_z = Zcoord(n)

        @views cumsum!(ll, d[index])
        ll .-= ll[ns+1]
        window .= exp.(.-(ll ./ (l[end] / length(l)) ./ (2ns + 1.0) .* 5.0) .^ 2)
        window ./= sum(window)
        unit_vector = sqrt((new_r - old_r)^2 + (new_z - old_z)^2)

        @inbounds for (k, i) in enumerate(index)
            nflux_r[i] += @. (new_r - old_r) / unit_vector * window[k] * W_per_trace / wall_s[i]
            nflux_z[i] += @. (new_z - old_z) / unit_vector * window[k] * W_per_trace / wall_s[i]
        end
    end

    # IMAS assignments
    dd = actor.dd
    dd.neutronics.first_wall.r = rwall
    dd.neutronics.first_wall.z = zwall
    ntt = resize!(dd.neutronics.time_slice)
    ntt.wall_loading.flux_r = nflux_r
    ntt.wall_loading.flux_z = nflux_z
    ntt.wall_loading.power = sqrt.(nflux_r .^ 2.0 .+ nflux_z .^ 2.0) .* wall_s

    # renormalize to ensure perfect power match
    norm = IMAS.fusion_neutron_power(dd.core_profiles.profiles_1d[]) / sum(ntt.wall_loading.power)
    ntt.wall_loading.flux_r .*= norm
    ntt.wall_loading.flux_z .*= norm
    ntt.wall_loading.power .*= norm

    # plot neutron wall loading cx
    if do_plot
        plot!(p, ntt.wall_loading; xlim=[minimum(rwall) * 0.9, maximum(rwall) * 1.1], title="First wall neutron loading")
        plot!(p, rwall, zwall; color=:black, label="", xlim=[minimum(rwall) * 0.9, maximum(rwall) * 1.1], title="")
        display(p)
    end

    return actor
end

mutable struct neutron_particle{T<:Float64}
    x::T
    y::T
    z::T
    δvx::T
    δvy::T
    δvz::T
end

function Rcoord(n::neutron_particle)
    return sqrt(n.x^2 + n.y^2)
end

function Zcoord(n::neutron_particle)
    return n.z
end

function define_neutrons(actor::ActorNeutronics)
    return define_neutrons(actor.dd, actor.par.N)
end

function define_neutrons(dd::IMAS.dd, N::Int)
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]

    # in-plasma mask
    r = eqt.profiles_2d[1].grid.dim1
    z = eqt.profiles_2d[1].grid.dim2
    R = [rr for rr in r, zz in z] # 2D
    Z = [zz for rr in r, zz in z] # 2D
    rz_lcfs = collect(zip(eqt.boundary.outline.r, eqt.boundary.outline.z))
    mask = [IMAS.PolygonOps.inpolygon((rr, zz), rz_lcfs) for rr in r, zz in z]

    # 2D neutron source
    # NOTE: we add power of neutrons coming from (D+D→He3+n) as if they were from (D+D→He4+n)
    # D+T   -> He4 (3.518 MeV)  + n (14.072 MeV) + 17.59MeV
    # D+D   -> He3 (0.8175 MeV) + n (2.4525 MeV) + 3.27MeV
    # This is OK for calculating wall loading but blanket will need to distinguish between these two to account for different energy
    neutron_source_1d = IMAS.D_T_to_He4_heating(cp1d) .* 4.0 .+ IMAS.D_D_to_He3_heating(cp1d) .* 3.0
    tmp = eqt.profiles_2d[1].psi .* (mask .== 1) .+ eqt.profiles_2d[1].psi[end] .* .*(mask .== 0) # to avoid extrapolation
    neutron_source_2d = IMAS.interp1d(cp1d.grid.psi, neutron_source_1d).(tmp) # W/m^3
    neutron_source_2d .*= mask .== 1 # set things to zero outsize of lcfs
    neutron_source_2d .*= R .* (2π * (r[2] - r[1]) * (z[2] - z[1])) # W

    # cumulative distribution function
    CDF = cumsum(neutron_source_2d[:])
    W_per_trace = CDF[end] / N
    CDF .= (CDF .- CDF[1]) ./ (CDF[end] - CDF[1])
    ICDF = IMAS.interp1d(CDF, Float64.(1:length(CDF)), :linear)

    # neutron structures
    neutrons = Vector{neutron_particle{Float64}}(undef, N)
    for k in eachindex(neutrons)
        dk = Int(ceil(ICDF(rand())))
        ϕ = rand() * 2π
        θv = rand() * 2π
        ϕv = acos(rand() * 2.0 - 1.0)
        Rk = R[dk]
        Zk = Z[dk]
        yk, xk = Rk .* sincos(ϕ)
        zk = Zk
        δvyk, δvxk = (sin(ϕv)) .* sincos(θv)
        δvzk = cos(ϕv)
        neutrons[k] = neutron_particle(xk, yk, zk, δvxk, δvyk, δvzk)
    end

    return neutrons, W_per_trace
end

@recipe function plot_neutrons(neutrons::Vector{neutron_particle{T}}, eqt::IMAS.equilibrium__time_slice) where {T<:Real}
    N = length(neutrons)
    r = eqt.profiles_2d[1].grid.dim1
    z = eqt.profiles_2d[1].grid.dim2

    # normalization weight to bring bin value in the unitary range
    plot_norm = 0.5 / (N * (r[2] - r[1]) * (z[2] - z[1]) / eqt.profiles_1d.area[end])

    @series begin
        seriestype --> :histogram2d
        nbins := (LinRange(minimum(r), maximum(r), length(r) - 1), LinRange(minimum(z), maximum(z), length(z) - 1))
        aspect_ratio := :equal
        grid := false
        weights := zeros(N) .+ plot_norm
        Rcoord.(neutrons), Zcoord.(neutrons)
    end
end

function define_wall(actor::ActorNeutronics; step::Float64=0.1)
    # resample wall and make sure it's clockwise (for COCOS = 11)
    eqt = actor.dd.equilibrium.time_slice[]
    wall = IMAS.first_wall(actor.dd.wall)
    rwall, zwall = IMAS.resample_2d_path(wall.r, wall.z; step, method=:linear)
    R0 = eqt.global_quantities.magnetic_axis.r
    Z0 = eqt.global_quantities.magnetic_axis.z
    IMAS.reorder_flux_surface!(rwall, zwall, R0, Z0; force_close=true)
    return rwall, zwall
end

function toroidal_intersection(r1::Real, z1::Real, r2::Real, z2::Real, x::Real, y::Real, z::Real, vx::Real, vy::Real, vz::Real)
    return toroidal_intersection(r1, z1, r2, z2, x, y, z, vx, vy, vz, vx^2 + vy^2, vz^2)
end

function toroidal_intersection(r1::Real, z1::Real, r2::Real, z2::Real, x::Real, y::Real, z::Real, vx::Real, vy::Real, vz::Real, v2::Real, vz2::Real)
    m = (z2 - z1) / (r2 - r1)
    z0 = z1 - m * r1
    m2 = m^2

    a = 0.0
    b = 0.0
    c = 0.0
    if isfinite(m)
        z_z0 = z - z0
        a = m2 * v2 - vz2
        b = 2 * (m2 * (vx * x + vy * y) - z_z0 * vz)
        c = m2 * (x^2 + y^2) - z_z0^2
    else
        a = v2
        b = 2 * (vx * x + vy * y)
        c = x^2 + y^2 - r1^2
    end

    t1 = -Inf
    t2 = -Inf
    if a == 0
        b != 0 && (t2 = -c / b)
    else
        nb_2a = -b / (2a)
        b2a_ca = nb_2a^2 - c / a
        if b2a_ca >= 0
            t_sq = sqrt(b2a_ca)
            t1 = nb_2a - t_sq
            t2 = nb_2a + t_sq
        end
    end

    t = Inf
    if t1 > 0
        zi = vz * t1 + z
        if zi > z1 && zi < z2
            t = t1
        end
    end
    if !isfinite(t) && t2 > 0
        zi = vz * t2 + z
        if zi > z1 && zi < z2
            t = t2
        end
    end

    return t
end

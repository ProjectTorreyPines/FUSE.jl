#= =============== =#
#  ActorNeutronics  #
#= =============== =#
Base.@kwdef mutable struct FUSEparameters__ActorNeutronics{T<:Real} <: ParametersActorPlasma{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    N::Entry{Int} = Entry{Int}("-", "Number of particles"; default=100000)
    do_plot::Entry{Bool} = act_common_parameters(do_plot=false)
end

mutable struct ActorNeutronics{D,P} <: SingleAbstractActor{D,P}
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

    neutrons, W_per_trace, dr, dz = define_neutrons(actor)
    if do_plot
        plot!(p, neutrons, actor.dd.equilibrium.time_slice[])
    end

    rwall, zwall = define_wall(actor)
    nflux_r, nflux_z, wall_s = IMAS.find_flux(neutrons, W_per_trace, rwall, zwall, dr,dz)

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


function define_neutrons(actor::ActorNeutronics)
    return define_neutrons(actor.dd, actor.par.N)
end

function define_neutrons(dd::IMAS.dd, N::Int)
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]
    source_1d = IMAS.D_T_to_He4_heating(cp1d) .* 4.0 .+ IMAS.D_D_to_He3_heating(cp1d) .* 3.0
    psi = cp1d.grid.psi
    neutrons, W_per_trace, dr, dz = IMAS.define_particles(eqt, psi, source_1d, N)

    return neutrons, W_per_trace, dr, dz
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

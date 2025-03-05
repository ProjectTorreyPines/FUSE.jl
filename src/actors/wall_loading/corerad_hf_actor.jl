#= ===================== =#
#  ActorCoreRadHeatFlux  #
#= ===================== =#
Base.@kwdef mutable struct FUSEparameters__ActorCoreRadHeatFlux{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    N::Entry{Int} = Entry{Int}("-", "Number of launched photons"; default=100000)
    r::Entry{Vector{T}} = Entry{Vector{T}}("m", "Vector of r at outermidplane"; default=T[])
    q::Entry{Vector{T}} = Entry{Vector{T}}("W m^-2", "Vector of parallel power density at outer midplane"; default=T[])
    levels::Entry{Union{Int,Vector}} =
        Entry{Union{Int,Vector}}("-", "If Int it defines number of levels in SOL, if vector it corresponds to the psi levels to build SOL"; default=20)
    merge_wall::Entry{Bool} = Entry{Bool}("-", "Merge dd.wall in mesh for the heat flux "; default=true)
    step::Entry{T} = Entry{T}("m", " Step for discretization of the default wall mesh (dd.wall)"; default=0.1)
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct ActorCoreRadHeatFlux{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorCoreRadHeatFlux{P}
    wall_heat_flux::Union{Nothing,IMAS.WallHeatFlux}
end

function ActorCoreRadHeatFlux(dd::IMAS.dd{D}, par::FUSEparameters__ActorCoreRadHeatFlux{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorCoreRadHeatFlux)
    par = par(kw...)
    return ActorCoreRadHeatFlux(dd, par, nothing)
end

"""
    ActorCoreRadHeatFlux(dd::IMAS.dd, act::ParametersAllActors; kw...)

Computes the heat flux on the wall due to the core radiation
"""
function ActorCoreRadHeatFlux(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorCoreRadHeatFlux(dd, act.ActorCoreRadHeatFlux; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorCoreRadHeatFlux{D,P}) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    # Compute wall mesh 
    # (Rwall, Zwall)        wall mesh (with intersections of SOL) -  m
    # s                     curvilinear abscissa computed from (Rwall, Zwall), clockwise starting at OMP - m
    Rwall, Zwall, _, _, s, _, _, _ = IMAS.mesher_heat_flux(dd; par.r, par.q, par.merge_wall, par.levels, par.step)

    # Parameters for heat flux due to core radiarion
    total_rad_source1d = IMAS.total_radiation_sources(dd.core_sources, cp1d; time0=dd.global_time)
    psi = cp1d.grid.psi
    source_1d = -total_rad_source1d.electrons.energy # minus sign because loss for dd.core_sources
    Prad_core = -total_rad_source1d.electrons.power_inside[end]

    # Compute the heat flux due to the core radiation - q_core_rad - W/m2
    q_core_rad = IMAS.core_radiation_heat_flux(eqt, psi, source_1d, par.N, Rwall, Zwall, Prad_core)

    HF = IMAS.WallHeatFlux{D}(; r=Rwall, z=Zwall, q_core_rad, s)
    actor.wall_heat_flux = HF

    #plot
    if par.do_plot
        ll = @layout [a{0.6w,0.9h} b{0.4w}]
        p = plot(; layout=ll, size=(1500, 500))
        plot!(p, HF; which_plot=:oneD, q=:core_radiation, subplot=1)
        sol = IMAS.sol(dd; levels=1)
        photons, W_per_trace, dr, dz = IMAS.define_particles(eqt, psi, source_1d, par.N)
        plot!(p, photons, actor.dd.equilibrium.time_slice[]; colorbar_entry=false, subplot=2)
        plot!(p, sol; subplot=2, line_z=nothing, color=:black)
        plot!(p, HF; q=:core_radiation, plot_type=:path, subplot=2)
        display(p)
    end

    return actor
end

function _finalize(actor::ActorCoreRadHeatFlux)
    # Finalize: work in progress - populate dd
    return actor
end
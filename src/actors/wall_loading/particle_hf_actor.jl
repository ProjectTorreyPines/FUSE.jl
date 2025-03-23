#= ===================== =#
#  ActorParticleHeatFlux  #
#= ===================== =#

Base.@kwdef mutable struct FUSEparameters__ActorParticleHeatFlux{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    r::Entry{Vector{T}} = Entry{Vector{T}}("m", "Vector of r at outermidplane"; default=T[])
    q::Entry{Vector{T}} = Entry{Vector{T}}("W m^-2", "Vector of parallel power density at outer midplane"; default=T[])
    levels::Entry{Union{Int,Vector{T}}} =
        Entry{Union{Int,Vector{T}}}("-", "If Int it defines number of levels in SOL, if vector it corresponds to the psi levels to build SOL"; default=20)
    merge_wall::Entry{Bool} = Entry{Bool}("-", "Merge dd.wall in mesh for the heat flux "; default=true)
    step::Entry{T} = Entry{T}("m", " Step for discretization of the default wall mesh (dd.wall)"; default=0.1)
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct ActorParticleHeatFlux{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorParticleHeatFlux{P}}
    wall_heat_flux::Union{Nothing,IMAS.WallHeatFlux}
end

function ActorParticleHeatFlux(dd::IMAS.dd{D}, par::FUSEparameters__ActorParticleHeatFlux{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorParticleHeatFlux)
    par = OverrideParameters(par; kw...)
    return ActorParticleHeatFlux(dd, par, nothing)
end

"""
    ActorParticleHeatFlux(dd::IMAS.dd, act::ParametersAllActors; kw...)

Computes the heat flux on the wall due to the charged particles
"""
function ActorParticleHeatFlux(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorParticleHeatFlux(dd, act.ActorParticleHeatFlux; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorParticleHeatFlux{D,P}) where {D<:Real, P<:Real}
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]

    # Compute wall mesh 
    # (Rwall, Zwall)        wall mesh (with intersections of SOL) -  m
    # (rwall, zwall)        evenly spaced wall mesh - m
    # s                     curvilinear abscissa computed from (Rwall, Zwall), clockwise starting at OMP - m
    # SOL                   list of OpenFieldLines used to compute (Rwall, Zwall) 
    # (r,q)                 Hypothesis of power density decay at omp for definition of SOL
    Rwall, Zwall, rwall, zwall, s, SOL, r, q = IMAS.mesher_heat_flux(dd; par.r, par.q, par.merge_wall, par.levels, par.step)

    # Compute the heat flux due to the influx of charged particles
    # q_part                 Heat flux due to particles perpendicular to the wall       - W/m2
    # q_parallel             Heat flux due to particles parallel to the magnetic field  - W/m2
    q_part, q_parallel = IMAS.particle_heat_flux(eqt, SOL, rwall, zwall, r, q; par.merge_wall)

    HF = IMAS.WallHeatFlux{D}(; r=Rwall, z=Zwall, q_part, q_parallel, s)
    actor.wall_heat_flux = HF

    #plot
    if par.do_plot
        ll = @layout [a{0.6w,0.9h} b{0.4w}]
        p = plot(; layout=ll, size=(1500, 500))
        plot!(p, HF; which_plot=:oneD, q=:particle, subplot=1)
        sol = IMAS.sol(dd; levels=100)
        plot!(p, sol; subplot=2, colorbar_entry=false, color=:grays)
        plot!(p, HF; q=:particle, plot_type=:path, subplot=2)
        display(p)
    end

    return actor
end

# Populate dd
function _finalize(actor::ActorParticleHeatFlux)
    # Finalize: work in progress - populate dd
    return actor
end
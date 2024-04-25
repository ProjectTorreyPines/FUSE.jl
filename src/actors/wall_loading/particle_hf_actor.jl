#= ===================== =#
#  ActorParticleHeatFlux  #
#= ===================== =#

Base.@kwdef mutable struct FUSEparameters__ActorParticleHeatFlux{T<:Real} <: ParametersActorPlasma{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    r::Entry{Vector{T}} = Entry{Vector{T}}("m", "Vector of r at outermidplane"; default=T[])
    q::Entry{Vector{T}} = Entry{Vector{T}}("W m^-2", "Vector of parallel power density at outer midplane"; default=T[])
    levels::Entry{Union{Int,Vector{T}}} = Entry{Union{Int,Vector{T}}}("-", "If Int it defines number of levels in SOL, if vector it corresponds to the psi levels to build SOL"; default=20)
    merge_wall::Entry{Bool} = Entry{Bool}("-", "Merge dd.wall in mesh for the heat flux "; default=true)
    step::Entry{T} = Entry{T}("m", " Step for discretization of the default wall mesh (dd.wall)"; default=0.1)
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct ActorParticleHeatFlux{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorParticleHeatFlux{P}
    function ActorParticleHeatFlux(dd::IMAS.dd{D}, par::FUSEparameters__ActorParticleHeatFlux{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorParticleHeatFlux)
        par = par(kw...)
        return new{D,P}(dd, par)
    end
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

function _step(actor::ActorParticleHeatFlux)
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]

    # Compute wall mesh 
    # (Rwall, Zwall)        wall mesh (with intersections of SOL) -  m
    # (rwall, zwall)        evenly spaced wall mesh - m
    # s                     curvilinear abscissa computed from (Rwall, Zwall), clockwise starting at OMP - m
    # SOL                   list of OpenFieldLines used to compute (Rwall, Zwall) 
    # (r,q)                 Hypothesis of power density decay at omp for definition of SOL
    Rwall, Zwall, rwall, zwall, s, SOL, r, q = IMAS.mesher_HF(dd; par.r, par.q, par.merge_wall, par.levels, par.step)

    # Compute the heat flux due to the influx of charged particles
    # q_part                 Heat flux due to particles perpendicular to the wall       - W/m2
    # q_parallel             Heat flux due to particles parallel to the magnetic field  - W/m2
    q_part, q_parallel = IMAS.particle_HF(eqt, SOL, rwall, zwall, r, q; par.merge_wall)

    HF = IMAS.WallHeatFlux(; r=Rwall, z=Zwall, q_part, q_parallel, s)

    #plot
    if par.do_plot
        font = 15
        _, psi_separatrix = IMAS.find_psi_boundary(dd)
        surface, _ = IMAS.flux_surface(eqt, psi_separatrix, :open)
        ll = @layout [a{0.6w,0.9h} b{0.4w}]
        p = plot(; layout=ll, size=(1500, 500))
        plot!(HF; which_plot=:oneD, q=:part, subplot=1, xtickfont= font,  ytickfont= font,guidefont= font, legendfont= font, titlefont = font+5)
        plot!(HF; q=:part, plot_type=:path, subplot=2, xtickfont= font,  ytickfont= font,guidefont= font, legendfont= font, colorbartitlefont = font)
        for surf in surface
            plot!(surf; color=:grey, subplot=2)
        end
        display(p)
    end

    return actor
end

# Populate dd
function _finalize(actor::ActorParticleHeatFlux)
    # Finalize: work in progress - populate dd
    return actor
end
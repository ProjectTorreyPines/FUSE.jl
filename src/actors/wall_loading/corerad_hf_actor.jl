#= ===================== =#
#  ActorCoreRadHeatFlux  #
#= ===================== =#

Base.@kwdef mutable struct FUSEparameters__ActorCoreRadHeatFlux{T<:Real} <: ParametersActorPlasma{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    N::Entry{Int} = Entry{Int}("-", "Number of launched photons"; default=500000)
    r::Entry{Vector{T}} = Entry{Vector{T}}("m", "Vector of r at outermidplane"; default=T[])
    q::Entry{Vector{T}} = Entry{Vector{T}}("W m^-2", "Vector of parallel power density at outer midplane"; default=T[])
    levels::Entry{Union{Int,Vector}} = Entry{Union{Int,Vector}}("-", "If Int it defines number of levels in SOL, if vector it corresponds to the psi levels to build SOL"; default=20)
    merge_wall::Entry{Bool} = Entry{Bool}("-", "Merge dd.wall in mesh for the heat flux "; default=true)
    step::Entry{T} = Entry{T}("m", " Step for discretization of the default wall mesh (dd.wall)"; default=0.1)
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct ActorCoreRadHeatFlux{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorCoreRadHeatFlux{P}
    function ActorCoreRadHeatFlux(dd::IMAS.dd{D}, par::FUSEparameters__ActorCoreRadHeatFlux{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorCoreRadHeatFlux)
        par = par(kw...)
        return new{D,P}(dd, par)
    end
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

function _step(actor::ActorCoreRadHeatFlux)
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    # Compute wall mesh 
    # (Rwall, Zwall)        wall mesh (with intersections of SOL) -  m
    # s                     curvilinear abscissa computed from (Rwall, Zwall), clockwise starting at OMP - m
    Rwall, Zwall, _, _, s, _, _, _ = IMAS.mesher_HF(dd; par.r, par.q, par.merge_wall, par.levels, par.step)

    # Parameters for heat flux due to core radiarion
    total_rad_source1d = IMAS.total_radiation_sources(dd.core_sources, cp1d)
    psi = dd.core_sources.source[1].profiles_1d[1].grid.psi
    source_1d = -total_rad_source1d.electrons.energy # minus sign because loss for dd.core_sources
    Prad_core = -total_rad_source1d.electrons.power_inside[end]

    # Compute the heat flux due to the core radiation - q_core_rad - W/m2
    q_core_rad = IMAS.core_radiation_HF(eqt, psi, source_1d, par.N, Rwall, Zwall, Prad_core)

    HF = IMAS.WallHeatFlux(; r=Rwall, z=Zwall, q_core_rad, s)

    #plot
    if par.do_plot
        font = 15
        _, psi_separatrix = IMAS.find_psi_boundary(dd)
        surface, _ = IMAS.flux_surface(eqt, psi_separatrix, :open)
        ll = @layout [a{0.6w,0.9h} b{0.4w}]
        p = plot(; layout=ll, size=(1500, 500))
        plot!(HF; which_plot=:oneD, q=:corerad, subplot=1, xtickfont= font,  ytickfont= font,guidefont= font, legendfont= font, titlefont = font+5)
        plot!(HF; q=:corerad, plot_type=:path, subplot=2, xtickfont= font,  ytickfont= font,guidefont= font, legendfont= font, colorbartitlefont = font)
        for surf in surface
            plot!(surf; color=:grey, subplot=2)
        end
        display(p)
    end
    return actor
end

function _finalize(actor::ActorCoreRadHeatFlux)
    # Finalize: work in progress - populate dd
    return actor
end
#= ===================== =#
#  ActorCoreRadHeatFlux  #
#= ===================== =#

# Definiton of act parameters
Base.@kwdef mutable struct FUSEparameters__ActorCoreRadHeatFlux{T<:Real} <: ParametersActorPlasma{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    N::Entry{Int} = Entry{Int}("-", "Number of launched photons"; default=500000)
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    r::Entry{AbstractVector{Float64}} = Entry{AbstractVector{Float64}}("m", "Vector of r at outermidplane"; default = Float64[])
    q::Entry{AbstractVector{Float64}} = Entry{AbstractVector{Float64}}("W m^-2", "Vector of parallel power density  at outer midplane"; default = Float64[])
    merge_wall::Entry{Bool} = Entry{Bool}("-","Merge dd.wall in mesh for the heat flux "; default = true)
    levels::Entry{Union{Int,AbstractVector}} = Entry{Union{Int,AbstractVector}}("-", "If Int it defines number of levels in SOL, if vector it corresponds to the psi levels to build SOL"; default = 20)
    step::Entry{Real} = Entry{Real}("m"," Step for discretization of the default wall mesh (dd.wall)"; default = 0.1)
end
# Definition of Actor struct
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
ActorCoreRadHeatFlux

Computes the heat flux on the wall due to the core radiation

"""

function ActorCoreRadHeatFlux(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorCoreRadHeatFlux(dd, act.ActorCoreRadHeatFlux; kw...)
    step(actor)
    finalize(actor)
    return actor
end

# Operate actor
function _step(actor::ActorCoreRadHeatFlux)
    dd = actor.dd
    par = actor.par

    do_plot::Bool = par.do_plot
    r::AbstractVector{Float64} = par.r
    q::AbstractVector{Float64} = par.q
    merge_wall::Bool = par.merge_wall
    step::Real = par.step
    levels::Union{Int,AbstractVector} = par.levels
    N::Int = par.N
   
    eqt = dd.equilibrium.time_slice[]

    # Compute wall mesh 
    # (Rwall, Zwall)        wall mesh (with intersections of SOL) -  m
    # s                     curvilinear abscissa computed from (Rwall, Zwall), clockwise starting at OMP - m
    Rwall, Zwall, _, _, s, _, _, _= IMAS.mesher_HF(dd, r = r, q = q, merge_wall = merge_wall, levels = levels, step = step)

    #Parameters for heat flux due to core radiarion
    psi       = dd.core_sources.source[1].profiles_1d[1].grid.psi
    source_1d = -IMAS.total_radiation_sources(dd).electrons.energy # minus sign because loss for dd.core_sources
    Prad_core = -IMAS.total_radiation_sources(dd).electrons.power_inside[end]

    # Compute the heat flux due to the core radiation - Qrad - W/m2
    Qrad = IMAS.core_radiation_HF(eqt, psi, source_1d, N, Rwall, Zwall, Prad_core)

    HF = IMAS.WallHeatFlux(r = Rwall, z = Zwall, q_core_rad = Qrad, s = s)

    #plot
    if do_plot
        _, psi_separatrix = IMAS.find_psi_boundary(dd)
        surface, _ = IMAS.flux_surface(eqt,psi_separatrix, :open)
        ll = @layout [a{0.6w,0.9h} b{0.4w}]
        p = plot(; layout=ll, size=(1000, 500))
        plot!(HF, which_plot = :oneD, q = :corerad; subplot =1)
        plot!(HF, q = :corerad, plot_type =  :path; subplot =2)
        for surf in surface
            plot!(surf, color = :grey; subplot =2)
        end
        display(plot!())
    end
    return actor
end

# Populate dd
function _finalize(actor::ActorCoreRadHeatFlux)
    # Finalize: work in progress - populate dd
    return actor
end
#= ===================== =#
#  ActorParticleHeatFlux  #
#= ===================== =#

# Definiton of act parameters
Base.@kwdef mutable struct FUSEparameters__ActorParticleHeatFlux{T<:Real} <: ParametersActorPlasma{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    r::Entry{AbstractVector{Float64}} = Entry{AbstractVector{Float64}}("m", "Vector of r at outermidplane"; default = Float64[])
    q::Entry{AbstractVector{Float64}} = Entry{AbstractVector{Float64}}("W m^-2", "Vector of parallel power density  at outer midplane"; default = Float64[])
    merge_wall::Entry{Bool} = Entry{Bool}("-","Merge dd.wall in mesh for the heat flux "; default = true)
    levels::Entry{Union{Int,AbstractVector}} = Entry{Union{Int,AbstractVector}}("-", "If Int it defines number of levels in SOL, if vector it corresponds to the psi levels to build SOL"; default = 20)
    step::Entry{Real} = Entry{Real}("m"," Step for discretization of the default wall mesh (dd.wall)"; default = 0.1)
end
# Definition of Actor struct
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
    ActorParticleHeatFlux

Computes the heat flux on the wall due to the charged particles

"""

function ActorParticleHeatFlux(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorParticleHeatFlux(dd, act.ActorParticleHeatFlux; kw...)
    step(actor)
    finalize(actor)
    return actor
end

# Operate actor
function _step(actor::ActorParticleHeatFlux)
    dd = actor.dd
    par = actor.par

    do_plot::Bool = par.do_plot
    r::AbstractVector{Float64} = par.r
    q::AbstractVector{Float64} = par.q
    merge_wall::Bool = par.merge_wall
    step::Real = par.step
    levels::Union{Int,AbstractVector} = par.levels
   
    eqt = dd.equilibrium.time_slice[]

    # Compute wall mesh 
    # (Rwall, Zwall)        wall mesh (with intersections of SOL) -  m
    # (rwall, zwall)        evenly spaced wall mesh - m
    # s                     curvilinear abscissa computed from (Rwall, Zwall), clockwise starting at OMP - m
    # SOL                   list of OpenFieldLines used to compute (Rwall, Zwall) 
    # (r,q)                 Hypothesis of power density decay at omp for definition of SOL
    Rwall, Zwall, rwall, zwall, s, SOL, r, q= IMAS.mesher_HF(dd, r = r, q = q, merge_wall = merge_wall, levels = levels, step = step)

    # Compute the heat flux due to the influx of charged particles
    # Qpart                 Heat flux due to particles perpendicular to the wall       - W/m2
    # Qpara                 Heat flux due to particles parallel to the magnetic field  - W/m2
    Qpart, Qpara = IMAS.particle_HF(eqt, SOL, rwall, zwall, r, q; merge_wall)


    HF = IMAS.WallHeatFlux(r = Rwall, z = Zwall, q_part = Qpart, q_parallel = Qpara, s = s)

    #plot
    if do_plot
        _, psi_separatrix = IMAS.find_psi_boundary(dd)
        surface, _ = IMAS.flux_surface(eqt,psi_separatrix, :open)
        ll = @layout [a{0.6w,0.9h} b{0.4w}]
        p = plot(; layout=ll, size=(1000, 500))
        plot!(HF, which_plot = :oneD, q = :part; subplot =1)
        plot!(HF, q = :part, plot_type =  :path; subplot =2)
        for surf in surface
            plot!(surf, color = :grey; subplot =2)
        end
        display(plot!())
    end
    return actor
end

# Populate dd
function _finalize(actor::ActorParticleHeatFlux)
    # Finalize: work in progress - populate dd
    return actor
end
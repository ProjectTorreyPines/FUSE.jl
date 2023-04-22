import PlasmaFacingSurfaces

#= ========================= =#
#  ActorPlasmaFacingSurfaces  #
#= ========================= =#

Base.@kwdef mutable struct FUSEparameters__ActorPlasmaFacingSurfaces{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    divertors::PlasmaFacingSurfaces.DivertorsDesignParameters{T} = PlasmaFacingSurfaces.DivertorsDesignParameters{T}()
    main_chamber_walls::PlasmaFacingSurfaces.MainChamberWallsDesignParameters{T} = PlasmaFacingSurfaces.MainChamberWallsDesignParameters{T}()
end

mutable struct ActorPlasmaFacingSurfaces <: ReactorAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorPlasmaFacingSurfaces
    pfs::Union{Nothing,PlasmaFacingSurfaces.PFSDesign}
end

"""
    ActorPlasmaFacingSurfaces(dd::IMAS.dd, act::ParametersAllActors; kw...)

Takes equilibrium and builds a first wall (ie. all plasma facing surfaces, including divertor)

Creates `main_chamber_wall` and `upper_divertor` and/or `lower_divertor` 2D descriptions to the `wall` IDS

A divertor is composed of three sub-systems:
- targets that are intersected by a strike point and that will receive the highest possible heat flux
- baffles that are mostly receiving heat fluxes from radiations and are largely parallel to the plasma parallel flow
- a dome that is exposed to the private flux region and that is connecting the inner and outer components
"""
function ActorPlasmaFacingSurfaces(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorPlasmaFacingSurfaces
    actor = ActorPlasmaFacingSurfaces(dd, par; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorPlasmaFacingSurfaces(dd::IMAS.dd, par::FUSEparameters__ActorPlasmaFacingSurfaces; kw...)
    logging_actor_init(ActorPlasmaFacingSurfaces)
    par = par(kw...)
    return ActorPlasmaFacingSurfaces(dd, par, nothing)
end

function _step(actor::ActorPlasmaFacingSurfaces)
    par = actor.par
    actor.pfs = PlasmaFacingSurfaces.PFSDesign(actor.dd.equilibrium.time_slice[])
    pfs_pars = PlasmaFacingSurfaces.PFSDesignParameters(WeakRef(nothing), :not_set, par.divertors, par.main_chamber_walls)
    actor.pfs(pfs_pars)
    return actor
end

function _finalize(actor::ActorPlasmaFacingSurfaces; kw...)
    PlasmaFacingSurfaces.export2ddwall(actor.dd, actor.pfs; kw...)
    PlasmaFacingSurfaces.export2ddbuild(actor.dd, actor.pfs; kw...)
    return actor
end
#= ========================= =#
#  ActorPlasmaFacingSurfaces #
#= ========================= =#

import PlasmaFacingSurfaces

#==========#

mutable struct ActorPlasmaFacingSurfaces <: PlasmaAbstractActor
    dd::IMAS.dd
    par::PlasmaFacingSurfaces.PFSDesignParameters
    pfs::Union{Nothing,PlasmaFacingSurfaces.PFSDesign}
end

"""
    ActorPlasmaFacingSurfaces(dd::IMAS.dd, act::ParametersAllActors; kw...)
Takes equilibrium and builds a first wall (ie. all plasma facing surfaces, including divertor)
Creates `main_chamber_wall` and `upper_divertor` and/or `lower_divertor`  2D descriptions to the `wall` IDS
A divertor is composed of three sub-systems:
- targets that are intersected by a strike point and that will receive the highest possible heat flux
- baffles that are mostly receiving heat fluxes from radiations and are largely parallel to the plasma parallel flow
- a dome that is exposed to the private flux region and that is connecting the inner and outer components
"""
function ActorPlasmaFacingSurfaces(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorPlasmaFacingSurfaces(kw...)
    actor = ActorPlasmaFacingSurfaces(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorPlasmaFacingSurfaces(dd::IMAS.dd, par::PlasmaFacingSurfaces.FUSEparameters__ActorPlasmaFacingSurfaces; kw...)
    logging_actor_init(ActorPlasmaFacingSurfaces)
    par = par(kw...)
    return ActorPlasmaFacingSurfaces(dd, par, nothing)
end

function _step(actor::ActorPlasmaFacingSurfaces)
    actor.pfs = PlasmaFacingSurfaces.PFSDesign(actor.dd.equilibrium.time_slice[])
    actor.pfs(actor.par)
    return actor
end

function _finalize(actor::ActorPlasmaFacingSurfaces; kw...)
    PlasmaFacingSurfaces.export2ddwall(actor.dd, actor.pfs; kw...)
    PlasmaFacingSurfaces.export2ddbuild(actor.dd, actor.pfs; kw...)
    actor
end
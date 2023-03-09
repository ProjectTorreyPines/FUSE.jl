#= ========================= =#
#  ActorPlasmaFacingSurfaces #
#= ========================= =#

import PlasmaFacingSurfaces

# const pfc_materials = [:tungsten, :SiC, :graphite]

Base.@kwdef mutable struct FUSEparameters__TargetBoundary{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    λ_plasma::Entry{T} = Entry(T, "m", "Plasma width"; default=0.003)
    β_plasma::Entry{T} = Entry(T, "-", "Plasma width decay factor"; default=5.0)
    Λ_buffer::Entry{T} = Entry(T, "m", "Buffer between plasma and wall"; default=0.05)
end

Base.@kwdef mutable struct FUSEparameters__DivertorTarget{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    #material::Switch{Symbol} = Switch(Symbol, pfc_materials, "-", "Target surface material"; default=:tungsten)
    cfr::FUSEparameters__TargetBoundary{T} = FUSEparameters__TargetBoundary{T}() # cfr target boundary point 
    pfr::FUSEparameters__TargetBoundary{T} = FUSEparameters__TargetBoundary{T}() # pfr target boundary point
    Λ_peak::Entry{T} = Entry(T, "m", "Position of the plasma peak with respect to separatrix "; default=0.0) # 
    type::Switch{Symbol} = Switch(Symbol, [:flat], "-", "Type of target"; default=:flat)
    θ_target::Entry{T} = Entry(T, "-", "Inclination of the target with respect to the separatrix normal"; default=60.0)
end

Base.@kwdef mutable struct FUSEparameters__CFRDivertorBaffle{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    #material::Switch{Symbol} = Switch(Symbol, pfc_materials, "-", "Baffle surface material"; default=:tungsten)
    method::Switch{Symbol} = Switch(Symbol, [:xpt_equator, :manual], "-", "Method to define the baffle design point"; default=:xpt_equator)
    𝓁_baffle::Entry{T} = Entry(T, "m", "Poloidal distance from target boundary point to baffle design point"; default=0.0)
    r_baffle::Entry{T} = Entry(T, "m", "Radial distance from separatrix to baffle design point"; default=0.0)
end

Base.@kwdef mutable struct FUSEparameters__PFRDivertorBaffle{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    #material::Switch{Symbol} = Switch(Symbol, pfc_materials, "-", "Baffle surface material"; default=:tungsten)
    method::Switch{Symbol} = Switch(Symbol, [:none, :manual], "-", "Method to define the baffle design point"; default=:none)
    𝓁_baffle::Entry{T} = Entry(T, "m", "Poloidal distance from target boundary point to baffle design point"; default=0.0)
    r_baffle::Entry{T} = Entry(T, "m", "Radial distance from separatrix to baffle design point"; default=0.0)
end

Base.@kwdef mutable struct FUSEparameters__DivertorLegParameters{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    target::FUSEparameters__DivertorTarget{T} = FUSEparameters__DivertorTarget{T}()
    baffle_cfr::FUSEparameters__CFRDivertorBaffle{T} = FUSEparameters__CFRDivertorBaffle{T}() # common flux region
    baffle_pfr::FUSEparameters__PFRDivertorBaffle{T} = FUSEparameters__PFRDivertorBaffle{T}() # private flux region
    l_leg::Entry{T} = Entry(T, "m", "Length of the divertor leg"; default=1.0)
    l_buffer::Entry{T} = Entry(T, "m", "Distance between the strike point and the closest vessel structure (e.g. HT shield)"; default=0.25)
    method::Switch{Symbol} = Switch(Symbol, [:manual, :buffer], "-", "Method to define the divertor leg length"; default=:manual)
end

Base.@kwdef mutable struct FUSEparameters__DivertorDomeParameters{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    #material::Switch{Symbol} = Switch(Symbol, pfc_materials, "-", "Dome surface material"; default=:tungsten)
    #𝓁_xpt::Entry{T} = Entry(T, "-", "Distance between top of the dome and xpt"; default=0.0)
    #α_xpt::Entry{T} = Entry(T, "-", "Fraction of distance between top of the dome and xpt"; default=0.0)
end

Base.@kwdef mutable struct FUSEparameters__DivertorParameters{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    inner::FUSEparameters__DivertorLegParameters{T} = FUSEparameters__DivertorLegParameters{T}()
    outer::FUSEparameters__DivertorLegParameters{T} = FUSEparameters__DivertorLegParameters{T}()
    dome::FUSEparameters__DivertorDomeParameters{T} = FUSEparameters__DivertorDomeParameters{T}()
end

Base.@kwdef mutable struct FUSEparameters__MainChamberWallBoundary{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    λ_plasma::Entry{T} = Entry(T, "m", "Plasma width"; default=0.003)
    β_plasma::Entry{T} = Entry(T, "-", "Plasma width decay factor"; default=5.0)
    Λ_buffer::Entry{T} = Entry(T, "m", "Buffer between plasma and wall"; default=0.05)
end

Base.@kwdef mutable struct FUSEparameters__MainChamberWallDesign{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    #material::Switch{Symbol} = Switch(Symbol, pfc_materials, "-", "First wall surface material"; default=:tungsten)
    mid_plane::FUSEparameters__MainChamberWallBoundary{T} = FUSEparameters__MainChamberWallBoundary{T}()
    method::Switch{Symbol} = Switch(Symbol, [:plasma_width, :conformal], "-", "Method to draw main chamber wall"; default=:plasma_width)
end

Base.@kwdef mutable struct FUSEparameters__DivertorsDesign{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    lower::FUSEparameters__DivertorParameters{T} = FUSEparameters__DivertorParameters{T}() #lower divertor parameters
    upper::FUSEparameters__DivertorParameters{T} = FUSEparameters__DivertorParameters{T}() #upper divertor parameters
end

Base.@kwdef mutable struct FUSEparameters__MainChamberWallsDesign{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    outer::FUSEparameters__MainChamberWallDesign{T} = FUSEparameters__MainChamberWallDesign{T}()
    inner::FUSEparameters__MainChamberWallDesign{T} = FUSEparameters__MainChamberWallDesign{T}()
end

Base.@kwdef mutable struct FUSEparameters__ActorPlasmaFacingSurfaces{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    divertors::FUSEparameters__DivertorsDesign{T} = FUSEparameters__DivertorsDesign{T}()
    main_chamber_walls::FUSEparameters__MainChamberWallsDesign{T} = FUSEparameters__MainChamberWallsDesign{T}()
end

#==========#

mutable struct ActorPlasmaFacingSurfaces <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorPlasmaFacingSurfaces
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

function ActorPlasmaFacingSurfaces(dd::IMAS.dd, par::FUSEparameters__ActorPlasmaFacingSurfaces; kw...)
    logging_actor_init(ActorPlasmaFacingSurfaces)
    par = par(kw...)
    return ActorPlasmaFacingSurfaces(dd, par, nothing)
end

function _step(actor::ActorPlasmaFacingSurfaces)
    actor.pfs = PlasmaFacingSurfaces.PFSDesign(actor.dd.equilibrium.time_slice[])
    actor.pfs(actor.par)
    return actor
end

function _finalize(actor::ActorPlasmaFacingSurfaces)
    PlasmaFacingSurfaces.export2ddwall(actor.dd.wall.description_2d, actor.pfs)
end
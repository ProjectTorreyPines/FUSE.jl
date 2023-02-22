#= ========================= =#
#  ActorPlasmaFacingSurfaces #
#= ========================= =#

import PlasmaFacingSurfaces


import PlasmaFacingSurfaces

const pfc_materials = [:tungsten, :SiC, :graphite]

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
    material::Switch{:Symbol} = Switch(Symbol, pfc_materials, "-", "Target surface material"; default=:tungsten)
    cfr::FUSEparameters__TargetBoundary = FUSEparameters__TargetBoundary() # cfr target boundary point 
    pfr::FUSEparameters__TargetBoundary = FUSEparameters__TargetBoundary() # pfr target boundary point
    Λ_peak::Entry{T} = Entry(T, "m", "Position of the plasma peak with respect to separatrix "; default=0.0) # 
    type::Switch{Symbol} = Switch(Symbol, [:flat], "-", "Type of target"; default=:flat)
    θ_target::Entry{T} = Entry(T, "-", "Inclination of the target with respect to the separatrix normal"; default=60.0)
end

Base.@kwdef mutable struct FUSEparameters__CFRDivertorBaffle{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    material::Switch{:Symbol} = Switch(Symbol, pfc_materials, "-", "Baffle surface material"; default=:tungsten)
    method::Switch{Symbol} = Switch(Symbol, [:xpt_equator, :manual], "-", "Method to define the baffle design point"; default=:xpt_equator)
    𝓁_baffle::Entry{T} = Entry(T, "m", "Poloidal distance from target boundary point to baffle design point"; default=0.0)
    r_baffle::Entry{T} = Entry(T, "m", "Radial distance from separatrix to baffle design point"; default=0.0)
end

Base.@kwdef mutable struct FUSEparameters__PFRDivertorBaffle{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    material::Switch{:Symbol} = Switch(Symbol, pfc_materials, "-", "Baffle surface material"; default=:tungsten)
    method::Switch{Symbol} = Switch(Symbol, [:none, :manual], "-", "Method to define the baffle design point"; default=:none)
    𝓁_baffle::Entry{T} = Entry(T, "m", "Poloidal distance from target boundary point to baffle design point"; default=0.0)
    r_baffle::Entry{T} = Entry(T, "m", "Radial distance from separatrix to baffle design point"; default=0.0)
end

Base.@kwdef mutable struct FUSEparameters__DivertorLegParameters{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    target::FUSEparameters__DivertorTarget = FUSEparameters__DivertorTarget()
    baffle_cfr::FUSEparameters__CFRDivertorBaffle = CFRDivertorBaffle() # common flux region
    baffle_pfr::FUSEparameters__PFRDivertorBaffle = PFRDivertorBaffle() # private flux region
    l_leg::Entry{T} = Entry(T, "m", "Length of the divertor leg"; default=1.0)
    d_heat_shield::Entry{T} = Entry(T, "m", "Length of the divertor leg"; default=0.25)
    method::Switch{:Symbol} = Switch(Symbol, [:manual, :d_heat_shield], "-", "Method to define the divertor leg length"; default=:manual)    
end

Base.@kwdef mutable struct FUSEparameters__DivertorDomeParameters{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    material::Switch{:Symbol} = Switch(Symbol, pfc_materials, "-", "Dome surface material"; default=:tungsten)
    𝓁_xpt::Entry{T} = Entry(T, "-", "Distance between top of the dome and xpt"; default=0.0)
    α_xpt::Entry{T} = Entry(T, "-", "Fraction of distance between top of the dome and xpt"; default=0.0)
end

Base.@kwdef mutable struct FUSEparameters__DivertorParameters{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    inner::FUSEparameters__DivertorLegParameters = FUSEparameters__DivertorLegParameters()
    outer::FUSEparameters__DivertorLegParameters = FUSEparameters__DivertorLegParameters()
    dome::FUSEparameters__DivertorDomeParameters = FUSEparameters__DivertorDomeParameters()
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
    material::Switch{:Symbol} = Switch(Symbol, pfc_materials, "-", "First wall surface material"; default=:tungsten)
    mid_plane::FUSEparameters__MainChamberWallBoundary = FUSEparameters__MainChamberWallBoundary()
    method::Switch{Symbol} = Switch(Symbol, [:plasma_width, :conformal], "-", "Method to draw main chamber wall"; default=:plasma_width)
end

Base.@kwdef mutable struct FUSEparameters__DivertorsDesign{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    lower::FUSEparameters__DivertorParameters = FUSEparameters__DivertorParameters() #lower divertor parameters
    upper::FUSEparameters__DivertorParameters = FUSEparameters__DivertorParameters() #upper divertor parameters
end

Base.@kwdef mutable struct FUSEparameters__MainChamberWallsDesign{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    outer::FUSEparameters__MainChamberWallDesign = FUSEparameters__MainChamberWallDesign()
    inner::FUSEparameters__MainChamberWallDesign = FUSEparameters__MainChamberWallDesign()
end

Base.@kwdef mutable struct FUSEparameters__ActorPlasmaFacingSurfaces{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    divertors::FUSEparameters__DivertorsDesign = FUSEparameters__DivertorsDesign()
    main_chamber_walls::FUSEparameters__MainChamberWallsDesign = FUSEparameters__MainChamberWallsDesign()
end

mutable struct ActorPlasmaFacingSurfaces <: PlasmaAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorPlasmaFacingSurfaces
    pfs::Union{Nothing,PlasmaFacingSurfaces.PFSDesign}
end

"""
    ActorPlasmaFacingSurfaces(dd::IMAS.dd, act::ParametersAllActors; kw...)

Takes equilibrium and builds a first wall (ie. all plasma facing surfaces, including divertor)

Creates `upper_divertor`, `lower_divertor`, and `main_chamber_wall` 2D descriptions to the `wall` IDS

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

function _step(actor::ActorPlasmaFacingSurfaces; do_plot=true)
    dd = actor.dd
    actor.pfs = PlasmaFacingSurfaces.PFSDesign(dd.equilibrium.time_slice[])
    actor.pfs(actor.par; do_plot=do_plot)
    PlasmaFacingSurfaces.export2ddwall(dd, actor.pfs)
    return actor
end

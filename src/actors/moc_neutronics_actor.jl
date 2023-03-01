using PyCall

Base.@kwdef mutable struct FUSEparameters__ActorMOCNeutronics{T} <: ParametersActor where {T<:Float64}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
end

mutable struct ActorMOCNeutronics <: ReactorAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorMOCNeutronics
    TBR::Float64
    leakage::Vector{Float64}
    thickness::Float64
    function ActorMOCNeutronics(dd::IMAS.dd, par::FUSEparameters__ActorMOCNeutronics; kw...)
        logging_actor_init(ActorMOCNeutronics)
        par = par(kw...)
        return new(dd, par)
    end
end

"""
    ActorMOCNeutronics(dd::IMAS.dd, act::ParametersAllActors; kw...)
    Runs 1D cylindrical neutron tranpsort calculations based on the method of characteristics
"""
function ActorMOCNeutronics(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorMOCNeutronics(kw...)
    actor = ActorMOCNeutronics(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorMOCNeutronics)
    data_path = "/home/mclaughlink/dev/mgxs_generator/data/500_cycles/" # need to find a good place for the nuclear data
    data_filename= "combined_mgxs_lib.h5"
    layer_thicknesses = [layer.thickness for layer in actor.dd.build.layer if layer.fs==-1] * 100 # to cm
    layer_materials = [replace(layer.material, " " => "-") for layer in actor.dd.build.layer if layer.fs==-1] # only doing lfs right now.
    replace!(layer_materials, "Vacuum" => "Air-(dry,-near-sea-level)") # approximation, need to make data point for vacuum
    py"""
    from MOCNeutronics import cylinders_openmoc # need to figure out how to install the right python
    def get_tbr(layer_thicknesses, layer_materials, data_path, data_filename):
        return cylinders_openmoc.run_openmoc(layer_thicknesses, layer_materials, data_path, data_filename)
    """
    tbr, leakage, thickness = py"get_tbr"(layer_thicknesses, layer_materials, data_path, data_filename)
    thickness = thickness/100 # back to meters
    actor.TBR = tbr
    actor.leakage = leakage
    actor.thickness = thickness
    return actor
end

function MOCNeutronicsOptimization(ini, act)
    dd = FUSE.init(ini,act)
    actor = FUSE.ActorMOCNeutronics(dd, act)
    objective_TBR = actor.TBR
    objective_shielding = sum(actor.leakage)
    objective_compactness = sum(actor.thickness)
    return [objective_TBR, objective_shielding, objective_compactness]
end
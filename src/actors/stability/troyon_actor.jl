import TroyonBetaNN

#= ================= =#
#  ActorTroyonBetaNN  #
#= ================= =#
Base.@kwdef mutable struct FUSEparameters__ActorTroyonBetaNN{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    verbose::Entry{Bool} = act_common_parameters(; verbose=false)
end

mutable struct ActorTroyonBetaNN{D,P} <: AbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorTroyonBetaNN{P}
    TD::TroyonBetaNN.Troyon_Data
end

function ActorTroyonBetaNN(dd::IMAS.dd{D}, par::FUSEparameters__ActorTroyonBetaNN{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorTroyonBetaNN)
    par = par(kw...)
    TD = TroyonBetaNN.load_predefined_Troyon_NN_Models();
    return ActorTroyonBetaNN{D,P}(dd, par, TD)
end

"""
    ActorTroyonBetaNN(dd::IMAS.dd, act::ParametersAllActors; kw...)

This actor evaluates the low-n no-wall ideal MHD stability
"""
function ActorTroyonBetaNN(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorTroyonBetaNN(dd, act.ActorTroyonBetaNN; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorTroyonBetaNN)
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]
    TroyonBetaNN.calculate_Troyon_beta_limits_for_a_given_time_slice(actor.TD, eqt; silence=true, par.verbose);
    return actor
end

function _finalize(actor::ActorTroyonBetaNN)
    dd = actor.dd
    par = actor.par

    mhd = resize!(dd.mhd_linear.time_slice; wipe=false)

    for MLP in actor.TD.MLPs
        mode = resize!(mhd.toroidal_mode, "perturbation_type.name" => "Troyon no-wall", "n_tor"=>MLP.n)
        mode.perturbation_type.description = "Troyon Beta NN (MLP) limit, no wall"
        mode.stability_metric = MLP.βₙ_limit
    end

    return actor
end

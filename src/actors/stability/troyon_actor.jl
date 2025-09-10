import TroyonBetaNN

#= ================= =#
#  ActorTroyonBetaNN  #
#= ================= =#
@actor_parameters_struct ActorTroyonBetaNN{T} begin
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
    verbose::Entry{Bool} = act_common_parameters(; verbose=false)
end

mutable struct ActorTroyonBetaNN{D,P} <: AbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorTroyonBetaNN{P}}
    TD::TroyonBetaNN.Troyon_Data
end

function ActorTroyonBetaNN(dd::IMAS.dd{D}, par::FUSEparameters__ActorTroyonBetaNN{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorTroyonBetaNN)
    par = OverrideParameters(par; kw...)
    TD = TroyonBetaNN.load_predefined_Troyon_NN_Models()
    return ActorTroyonBetaNN{D,P}(dd, par, TD)
end

"""
    ActorTroyonBetaNN(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates ideal MHD stability limits using neural network models for low-n no-wall modes.
The actor calculates the Troyon beta limit (βN) for different toroidal mode numbers using 
machine learning models trained on MHD stability analysis.

The calculation is disabled for negative triangularity plasmas. Results are stored in 
`dd.mhd_linear.time_slice` with separate entries for each toroidal mode number.

!!! note

    Stores data in `dd.mhd_linear`
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
    if eqt.boundary.triangularity < 0.0
        @warn("ActorTroyonBetaNN is disabled for negative triangularity plasmas")
    else
        TroyonBetaNN.calculate_Troyon_beta_limits_for_a_given_time_slice(actor.TD, eqt; silence=true, par.verbose)
    end

    return actor
end

function _finalize(actor::ActorTroyonBetaNN)
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]
    if eqt.boundary.triangularity >= 0.0
        mhd = resize!(dd.mhd_linear.time_slice; wipe=false)
        for MLP in actor.TD.MLPs
            mode = resize!(mhd.toroidal_mode, "perturbation_type.name" => "Troyon no-wall", "n_tor" => MLP.n)
            mode.perturbation_type.description = "Troyon Beta NN (MLP) limit, no wall"
            mode.stability_metric = MLP.βₙ_limit
        end
    end

    return actor
end

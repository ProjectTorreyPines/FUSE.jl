import VacuumFields

#= ====================== =#
#  ActorVerticalStability  #
#= ====================== =#
Base.@kwdef mutable struct FUSEparameters__ActorVerticalStability{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    #== actor parameters ==#
    wall_precision::Entry{Float64} = Entry{Float64}("-", "Precision for making wall quadralaterals"; default=0.1)
    wall_max_seg_length::Entry{Float64} = Entry{Float64}("-", "Maximum segment length for making wall quadralaterals"; default=0.5)
    do_plot::Entry{Bool} = act_common_parameters(do_plot=false)
end

mutable struct ActorVerticalStability{D,P} <: ReactorAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorVerticalStability{P}
    stability_margin::D
    normalized_growth_rate::D
    function ActorVerticalStability(dd::IMAS.dd{D}, par::FUSEparameters__ActorVerticalStability{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorVerticalStability)
        par = par(kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorVerticalStability(dd::IMAS.dd, act::ParametersAllActors; kw...)

Compute vertical stability metrics
"""
function ActorVerticalStability(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorVerticalStability(kw...) # this makes a local copy of `act.ActorVerticalStability` and overrides it with keywords that the user may have passed
    actor = ActorVerticalStability(dd, par) # instantiate the actor (see function below)
    step(actor)                # run the actor
    finalize(actor)            # finalize
    return actor
end

"""
    _step(actor::ActorVerticalStability)

Compute vertical stability metrics
"""
function _step(actor::ActorVerticalStability)
    par = actor.par
    dd = actor.dd
    bd = dd.build
    eqt = dd.equilibrium.time_slice[]
    Ip = eqt.global_quantities.ip
    active_coils = IMAS_pf_active__coils(dd; green_model=:realistic)
    quads = wall_quads(bd, par.wall_precision, par.wall_max_seg_length)
    passive_coils =  [VacuumFields.QuadCoil(R, Z) for (R, Z) in quads]
    coils = vcat(active_coils, passive_coils)
    for coil in coils
        coil.resistance = 1.0
    end
    image = VacuumFields.Image(dd)

    actor.stability_margin = VacuumFields.stability_margin(image, coils, Ip)
    actor.normalized_growth_rate = 0.0 #VacuumFields.normalized_growth_rate(image, coils, Ip)
    return actor
end

"""
    _finalize(actor::ActorVerticalStability)

Store vertical stability metrics
"""
function _finalize(actor::ActorVerticalStability)
    dd = actor.dd
    mhd = resize!(dd.mhd_linear.time_slice)
    resize!(mhd.toroidal_mode, 2)

    # Stability margin
    mode = mhd.toroidal_mode[1]
    mode.perturbation_type.description = "Vertical stability margin, > 0.15 for stability (N.B., not in Hz)"
    mode.perturbation_type.name = "m_s"
    mode.n_tor = 0
    mode.growthrate = actor.stability_margin # not in Hz

    # Normalized growth rate
    mode = mhd.toroidal_mode[2]
    mode.perturbation_type.description = "Normalized vertical growth rate, < 10 for stability (N.B., not in Hz)"
    mode.perturbation_type.name = "γτ"
    mode.n_tor = 0
    mode.growthrate = actor.normalized_growth_rate # not in Hz

    return actor
end
import VacuumFields

#= ====================== =#
#  ActorVerticalStability  #
#= ====================== =#
Base.@kwdef mutable struct FUSEparameters__ActorVerticalStability{T<:Real} <: ParametersActorBuild{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    #== actor parameters ==#
    wall_precision::Entry{Float64} = Entry{Float64}("-", "Precision for making wall quadralaterals"; default=0.1)
    wall_max_seg_length::Entry{Float64} = Entry{Float64}("-", "Maximum segment length for making wall quadralaterals"; default=0.5)
    default_passive_material::Entry{Symbol} = Entry{Symbol}("-", "Default material to use for poorly defined vacuum vessel"; default=:steel)
    do_plot::Entry{Bool} = act_common_parameters(do_plot=false)
end

mutable struct ActorVerticalStability{D,P} <: SingleAbstractActor{D,P}
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
    active_coils = IMAS_pf_active__coils(dd; green_model=:quad)

    # Defaults
    actor.stability_margin = NaN
    actor.normalized_growth_rate = NaN

    if all(coil.current == 0.0 for coil in active_coils)
        @warn "Active coils have no current. Can't compute vertical stability metrics"
        return actor
    end

    # The vacuum vessel can have multiple layers
    # Find all the bounding ones and turn the area in between into quads
    kvessel = IMAS.get_build_indexes(bd.layer; type=_vessel_, fs=_lfs_)
    if isempty(kvessel)
        @warn "No vessel found. Can't compute vertical stability metrics"
        return actor
    end
    kout = kvessel[end]
    kin = kvessel[1] - 1
    quads = layer_quads(bd.layer[kin], bd.layer[kout], par.wall_precision, par.wall_max_seg_length)

    passive_coils =  [VacuumFields.QuadCoil(R, Z) for (R, Z) in quads]

    # Compute resistance based on material resistivity & geometry of coil
    # N.B.: this just takes the material from the outermost build layer;
    #       does not account for toroidal breaks, heterogeneous materials,
    #       or builds with "water" vacuum vessels
    mat_vv = Material(bd.layer[kout].material)
    if ismissing(mat_vv) || ismissing(mat_vv.electrical_conductivity)
        mat_vv = Material(par.default_passive_material)
    end
    eta = 1.0 / mat_vv.electrical_conductivity(temperature=0.0)
    if !ismissing(eta)
        for coil in passive_coils
            coil.resistance = VacuumFields.resistance(coil, eta)
        end
    end

    image = VacuumFields.Image(eqt)

    coils = vcat(active_coils, passive_coils)

    actor.stability_margin = VacuumFields.stability_margin(image, coils, Ip)

    for (k, coil) in enumerate(active_coils)
        if coil.resistance <= 0.0
            @warn "Active coil #$(k) has invalid resistance: $(coil.resistance). Can't compute normalized growth rate.\nOffending coil: $(repr(coil))"
            return actor
        end
    end
    for (k, coil) in enumerate(passive_coils)
        if coil.resistance <= 0.0
            @warn "Passive coil #$(k) has invalid resistance: $(coil.resistance). Can't compute normalized growth rate.\nOffending coil: $(repr(coil))"
            return actor
        end
    end

    _, _, actor.normalized_growth_rate = VacuumFields.normalized_growth_rate(image, coils, Ip)

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
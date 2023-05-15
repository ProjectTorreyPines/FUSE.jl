#= ========== =#
#  flux-swing #
#= ========== =#
Base.@kwdef mutable struct FUSEparameters__ActorFluxSwing{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    operate_oh_at_j_crit::Entry{Bool} = Entry(Bool, "-", """
If `true` it makes the OH operate at its current limit (within specified `j_tolerance`).
The flattop duration and maximum toroidal magnetic field follow from that.
Otherwise we evaluate what is the current needed for a given flattop duration,
which may or may not exceed the OH critical current limit.""";
        default=true
    )
    j_tolerance::Entry{T} = Entry(T, "-", "Tolerance on the OH current limit"; default=0.4)
end

mutable struct ActorFluxSwing <: ReactorAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorFluxSwing
    function ActorFluxSwing(dd::IMAS.dd, par::FUSEparameters__ActorFluxSwing; kw...)
        logging_actor_init(ActorFluxSwing)
        par = par(kw...)
        return new(dd, par)
    end
end

"""
    ActorFluxSwing(dd::IMAS.dd, act::ParametersAllActors; kw...)

Depending on `operate_oh_at_j_crit`
* true => Evaluate the OH current limits, and evaluate flattop duration from that.
* false => Evaluate what are the currents needed for a given flattop duration. This may or may not exceed the OH current limits.

OH flux consumption based on:
* rampup estimate based on Ejima coefficient
* flattop consumption
* vertical field from PF coils

!!! note
    Stores data in `dd.build.flux_swing`, `dd.build.tf`, and `dd.build.oh`
"""
function ActorFluxSwing(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorFluxSwing(dd, act.ActorFluxSwing; kw...)
    step(actor)
    finalize(actor)
    return actor
end

"""
    _step(actor::ActorFluxSwing)
"""
function _step(actor::ActorFluxSwing)
    bd = actor.dd.build
    par = actor.par

    requirements = actor.dd.requirements
    eq = actor.dd.equilibrium
    eqt = eq.time_slice[]
    cp = actor.dd.core_profiles
    cp1d = cp.profiles_1d[]

    # oh
    bd.flux_swing.rampup = rampup_flux_estimates(eqt, cp)
    bd.flux_swing.pf = pf_flux_estimates(eqt)
    if par.operate_oh_at_j_crit
        oh_maximum_J_B!(bd; par.j_tolerance)
        bd.flux_swing.flattop = flattop_flux_estimates(bd) # flattop flux based on available current
    else
        bd.flux_swing.flattop = flattop_flux_estimates(requirements, cp1d) # flattop flux based on requirements duration
        oh_required_J_B!(bd)
    end

    # flattop duration
    flattop_duration!(bd, cp1d)

    # tf
    tf_required_J_B!(bd, eq)

    return actor
end

"""
    rampup_flux_estimates(eqt::IMAS.equilibrium__time_slice, cp::IMAS.core_profiles)

Estimate OH flux requirement during rampup, where eqt is supposed to be the equilibrium right at the end of the rampup phase,
beginning of flattop and core_profiles is only used to get core_profiles.global_quantities.ejima
"""
function rampup_flux_estimates(eqt::IMAS.equilibrium__time_slice, cp::IMAS.core_profiles)
    ###### what equilibrium time-slice should we use to evaluate rampup flux requirements?

    # from IMAS dd to local variables
    majorRadius = eqt.boundary.geometric_axis.r
    minorRadius = eqt.boundary.minor_radius
    elongation = eqt.boundary.elongation
    plasmaCurrent = eqt.global_quantities.ip / 1E6 # in [MA]
    li = eqt.global_quantities.li_3 # what li
    ejima = @ddtime cp.global_quantities.ejima

    # evaluate plasma inductance
    plasmaInductanceInternal = 0.4 * 0.5 * pi * majorRadius * li
    plasmaInductanceExternal = 0.4 * pi * majorRadius * (log(8.0 * majorRadius / minorRadius / sqrt(elongation)) - 2.0)
    plasmaInductanceTotal = plasmaInductanceInternal + plasmaInductanceExternal

    # estimate rampup flux requirement
    rampUpFlux = (ejima * 0.4 * pi * majorRadius + plasmaInductanceTotal) * plasmaCurrent

    return abs(rampUpFlux)
end

"""
    flattop_flux_estimates(requirements::IMAS.requirements, cp1d::IMAS.core_profiles__profiles_1d)

Estimate OH flux requirement during flattop

NOTES:
* If j_ohmic profile is missing then steady state ohmic profile is assumed
"""
function flattop_flux_estimates(requirements::IMAS.requirements, cp1d::IMAS.core_profiles__profiles_1d)
    return abs(integrate(cp1d.grid.area, cp1d.j_ohmic ./ cp1d.conductivity_parallel)) * (requirements.flattop_duration) # V*s
end

"""
    flattop_flux_estimates(bd::IMAS.build; double_swing::Bool=true)

OH flux given its bd.oh.max_b_field and radial build
"""
function flattop_flux_estimates(bd::IMAS.build; double_swing::Bool=true)
    OH = IMAS.get_build(bd, type=_oh_)
    innerSolenoidRadius = OH.start_radius
    outerSolenoidRadius = OH.end_radius
    magneticFieldSolenoidBore = bd.oh.max_b_field
    RiRo_factor = innerSolenoidRadius / outerSolenoidRadius
    totalOhFluxReq = magneticFieldSolenoidBore / 3.0 * pi * outerSolenoidRadius^2 * (RiRo_factor^2 + RiRo_factor + 1.0) * (double_swing ? 2 : 1)
    bd.flux_swing.flattop = totalOhFluxReq - bd.flux_swing.rampup - bd.flux_swing.pf
end

"""
    pf_flux_estimates(eqt::IMAS.equilibrium__time_slice)

Estimate vertical field from PF coils and its contribution to flux swing,
where `eqt` is supposed to be the equilibrium right at the end of the rampup phase, beginning of flattop
"""
function pf_flux_estimates(eqt::IMAS.equilibrium__time_slice)
    # from IMAS dd to local variables
    majorRadius = eqt.boundary.geometric_axis.r
    minorRadius = eqt.boundary.minor_radius
    elongation = eqt.boundary.elongation
    plasmaCurrent = eqt.global_quantities.ip / 1E6 # in [MA]
    betaP = eqt.global_quantities.beta_pol
    li = eqt.global_quantities.li_3 # what li does Stambaugh FST 2011 use?

    # estimate vertical field and its contribution to flux swing
    verticalFieldAtCenter = 0.1 * plasmaCurrent / majorRadius * (log(8.0 * majorRadius / (minorRadius * sqrt(elongation))) - 1.5 + betaP + 0.5 * li)
    fluxFromVerticalField = 0.8 * verticalFieldAtCenter * pi * (majorRadius^2 - (majorRadius - minorRadius)^2)

    return -abs(fluxFromVerticalField)
end

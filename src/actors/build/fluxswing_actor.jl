#= ========== =#
#  flux-swing #
#= ========== =#
@actor_parameters_struct ActorFluxSwing{T} begin
    operate_oh_at_j_crit::Entry{Bool} = Entry{Bool}("-", """
If `true` it makes the OH operate at its current limit (within specified dd.requirements.coil_j_margin`).
The flattop duration and maximum toroidal magnetic field follow from that.
Otherwise we evaluate what is the current needed for dd.requirements.flattop_duration,
which may or may not exceed the OH critical current limit.
If dd.requirements.flattop_duration is not set, then operate_oh_at_j_crit is assumed.""";
        default=false
    )
end

mutable struct ActorFluxSwing{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorFluxSwing{P}}
    function ActorFluxSwing(dd::IMAS.dd{D}, par::FUSEparameters__ActorFluxSwing{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorFluxSwing)
        par = OverrideParameters(par; kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorFluxSwing(dd::IMAS.dd, act::ParametersAllActors; kw...)

Analyzes OH coil flux swing capability and determines operating limits based on current density constraints
and flattop duration requirements. Calculates flux consumption during different operational phases
and determines maximum achievable flattop duration or required current levels.

# Operating modes (controlled by `operate_oh_at_j_crit`)
- `true`: Operates OH at critical current limit, calculates achievable flattop duration
- `false`: Uses required flattop duration, calculates needed currents (may exceed limits)

# Flux consumption analysis
- **Rampup flux**: Estimated using Ejima coefficient and plasma inductance calculations
- **Flattop flux**: Based on resistive current drive requirements from Ohmic heating  
- **PF flux**: Vertical field contribution from poloidal field coils

# Physics models
- Plasma inductance (internal + external components) for flux estimation
- Resistive flux consumption using conductivity profiles and current density
- Vertical field requirements based on plasma geometry and pressure

# Key outputs
- OH flux swing components (`dd.build.flux_swing.rampup/flattop/pf`)
- OH coil current and field limits (`dd.build.oh.max_j/max_b_field`)
- TF coil current and field requirements (`dd.build.tf.max_j/max_b_field`)  
- Achievable flattop duration (`dd.build.oh.flattop_duration`)

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

Performs flux swing analysis by calculating OH flux consumption components:
- Calls `rampup_flux_estimates()` using Ejima coefficient and plasma inductance
- Calls `pf_flux_estimates()` for vertical field contribution  
- Either calculates flattop flux from required duration or operates at current limit
- Updates OH and TF coil current/field requirements based on analysis results
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
    if !ismissing(requirements, :flattop_duration) && !par.operate_oh_at_j_crit
        bd.flux_swing.flattop = flattop_flux_estimates(requirements, cp1d; requirements.coil_j_margin) # flattop flux based on requirements duration
        oh_required_J_B!(bd)
    else
        oh_maximum_J_B!(bd; requirements.coil_j_margin)
        bd.flux_swing.flattop = flattop_flux_estimates(bd) # flattop flux based on available current
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
    flattop_flux_estimates(requirements::IMAS.requirements, cp1d::IMAS.core_profiles__profiles_1d; coil_j_margin::Float64)

Estimate OH flux requirement during flattop
"""
function flattop_flux_estimates(requirements::IMAS.requirements, cp1d::IMAS.core_profiles__profiles_1d; coil_j_margin::Float64)
    j_ohmic = cp1d.j_ohmic
    conductivity_parallel = cp1d.conductivity_parallel
    f = (k, x) -> j_ohmic[k] / conductivity_parallel[k]
    return abs(trapz(cp1d.grid.area, f)) * requirements.flattop_duration * (1.0 + coil_j_margin) # V*s
end

"""
    flattop_flux_estimates(bd::IMAS.build; double_swing::Bool=true)

OH flux given its bd.oh.max_b_field and radial build
"""
function flattop_flux_estimates(bd::IMAS.build; double_swing::Bool=true)
    OH = IMAS.get_build_layer(bd.layer; type=_oh_)
    innerSolenoidRadius = OH.start_radius
    outerSolenoidRadius = OH.end_radius
    magneticFieldSolenoidBore = bd.oh.max_b_field
    RiRo_factor = innerSolenoidRadius / outerSolenoidRadius
    totalOhFluxReq = magneticFieldSolenoidBore / 3.0 * pi * outerSolenoidRadius^2 * (RiRo_factor^2 + RiRo_factor + 1.0) * (double_swing ? 2 : 1)
    return bd.flux_swing.flattop = totalOhFluxReq - bd.flux_swing.rampup - bd.flux_swing.pf
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

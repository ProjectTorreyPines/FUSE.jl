#= ========= =#
#  OH magnet  #
#= ========= =#
"""
    oh_maximum_J_B!(bd::IMAS.build; j_tolerance)

Evaluate maxium OH current density and magnetic field for given geometry and technology

NOTES:
* Equations from GASC (Stambaugh FST 2011)
* Also relevant: `Engineering design solutions of flux swing with structural requirements for ohmic heating solenoids` Smith, R. A. September 30, 1977
"""
function oh_maximum_J_B!(bd::IMAS.build; j_tolerance)
    OH = IMAS.get_build(bd, type=_oh_)
    innerSolenoidRadius = OH.start_radius
    outerSolenoidRadius = OH.end_radius

    # find maximum superconductor critical_j given self-field
    function max_J_OH(x)
        currentDensityOH = abs(x[1])
        magneticFieldSolenoidBore = currentDensityOH / 1E6 * (0.4 * pi * outerSolenoidRadius * (1.0 - innerSolenoidRadius / outerSolenoidRadius))
        critical_j = coil_J_B_crit(magneticFieldSolenoidBore, bd.oh.technology)[1]
        # do not use relative error here. Absolute error tells optimizer to lower currentDensityOH if critical_j==0
        return abs(critical_j - currentDensityOH * (1.0 + j_tolerance))
    end
    res = Optim.optimize(max_J_OH, 0.0, 1E9, Optim.GoldenSection(), rel_tol=1E-3)

    # solenoid maximum current and field
    bd.oh.max_j = abs(res.minimizer[1])
    bd.oh.max_b_field = bd.oh.max_j / 1E6 * (0.4 * pi * outerSolenoidRadius * (1.0 - innerSolenoidRadius / outerSolenoidRadius))
    bd.oh.critical_j, bd.oh.critical_b_field = coil_J_B_crit(bd.oh.max_b_field, bd.oh.technology)
end

"""
    oh_required_J_B!(bd::IMAS.build; double_swing::Bool=true)

Evaluate OH current density and B_field required for given rampup and flattop
NOTES:
* Equations from GASC (Stambaugh FST 2011)
* Also relevant: `Engineering design solutions of flux swing with structural requirements for ohmic heating solenoids` Smith, R. A. September 30, 1977
"""
function oh_required_J_B!(bd::IMAS.build; double_swing::Bool=true)
    OH = IMAS.get_build(bd, type=_oh_)
    innerSolenoidRadius = OH.start_radius
    outerSolenoidRadius = OH.end_radius

    totalOhFluxReq = bd.flux_swing.rampup + bd.flux_swing.flattop + bd.flux_swing.pf

    # Calculate magnetic field at solenoid bore required to match flux swing request
    RiRo_factor = innerSolenoidRadius / outerSolenoidRadius
    magneticFieldSolenoidBore = 3.0 * totalOhFluxReq / pi / outerSolenoidRadius^2 / (RiRo_factor^2 + RiRo_factor + 1.0) / (double_swing ? 2 : 1)
    currentDensityOH = magneticFieldSolenoidBore / (0.4 * pi * outerSolenoidRadius * (1 - innerSolenoidRadius / outerSolenoidRadius))

    # minimum requirements for OH
    bd.oh.max_b_field = magneticFieldSolenoidBore
    bd.oh.max_j = currentDensityOH * 1E6
    bd.oh.critical_j, bd.oh.critical_b_field = coil_J_B_crit(bd.oh.max_b_field, bd.oh.technology)
end

"""
    flattop_duration!(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d; double_swing::Bool=true)

Estimate OH flux requirement during flattop (if j_ohmic profile is missing then steady state ohmic profile is assumed)
"""
function flattop_duration!(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice, cp1d::IMAS.core_profiles__profiles_1d; double_swing::Bool=true)
    OH = IMAS.get_build(bd, type=_oh_)
    innerSolenoidRadius = OH.start_radius
    outerSolenoidRadius = OH.end_radius

    # estimate oh flattop flux and duration
    RiRo_factor = innerSolenoidRadius / outerSolenoidRadius
    totalOhFlux = bd.oh.max_b_field * (pi * outerSolenoidRadius^2 * (RiRo_factor^2 + RiRo_factor + 1.0) * (double_swing ? 2 : 1)) / 3.0
    bd.flux_swing.flattop = totalOhFlux - bd.flux_swing.rampup - bd.flux_swing.pf
    bd.oh.flattop_duration = bd.flux_swing.flattop / abs(integrate(cp1d.grid.area, cp1d.j_ohmic ./ cp1d.conductivity_parallel))
end

#= ========= =#
#  TF magnet  #
#= ========= =#
"""
    tf_maximum_J_B!(bd::IMAS.build; j_tolerance)

Evaluate maxium TF current density and magnetic field for given geometry and technology
"""
function tf_maximum_J_B!(bd::IMAS.build; j_tolerance)
    hfsTF = IMAS.get_build(bd, type=_tf_, fs=_hfs_)
    TF_cx_area = hfsTF.thickness * bd.tf.wedge_thickness

    # find maximum superconductor critical_j given self-field
    function max_J_TF(x)
        currentDensityTF = abs(x[1])
        current_TF = currentDensityTF * TF_cx_area
        max_b_field = current_TF / hfsTF.end_radius / 2pi * constants.μ_0 * bd.tf.coils_n
        critical_j = coil_J_B_crit(max_b_field, bd.tf.technology)[1]
        # do not use relative error here. Absolute error tells optimizer to lower currentDensityTF if critical_j==0
        return abs(critical_j - currentDensityTF * (1.0 + j_tolerance))
    end
    res = Optim.optimize(max_J_TF, 0.0, 1E9, Optim.GoldenSection(), rel_tol=1E-3)

    # tf maximum current and field
    bd.tf.max_j = abs(res.minimizer[1])
    current_TF = bd.tf.max_j * TF_cx_area
    bd.tf.max_b_field = current_TF / hfsTF.end_radius / 2pi * constants.μ_0 * bd.tf.coils_n
    bd.tf.critical_j, bd.tf.critical_b_field = coil_J_B_crit(bd.tf.max_b_field, bd.tf.technology)
end

"""
    tf_required_J_B!(bd::IMAS.build)

Evaluate TF current density given a B_field
"""
function tf_required_J_B!(bd::IMAS.build, eq::IMAS.equilibrium)
    hfsTF = IMAS.get_build(bd, type=_tf_, fs=_hfs_)
    lfsTF = IMAS.get_build(bd, type=_tf_, fs=_lfs_)
    B0 = abs(maximum(eq.vacuum_toroidal_field.b0))
    R0 = (hfsTF.end_radius + lfsTF.start_radius) / 2.0

    # current in the TF coils
    current_TF = B0 * R0 * 2pi / constants.μ_0 / bd.tf.coils_n
    TF_cx_area = hfsTF.thickness * bd.tf.wedge_thickness

    bd.tf.max_b_field = B0 * R0 / hfsTF.end_radius
    bd.tf.max_j = current_TF / TF_cx_area
    bd.tf.critical_j, bd.tf.critical_b_field = coil_J_B_crit(bd.tf.max_b_field, bd.tf.technology)
end

#= ========== =#
#  flux-swing #
#= ========== =#
mutable struct ActorFluxSwing <: ReactorAbstractActor
    dd::IMAS.dd
    par::ParametersActor
    operate_at_j_crit::Bool
    j_tolerance::Real
end

function ParametersActor(::Type{Val{:ActorFluxSwing}})
    par = ParametersActor(nothing)
    par.operate_at_j_crit = Entry(
        Bool,
        "",
        """
Makes the OH and TF operate at their current limit (within specified `j_tolerance`).
The flattop duration and maximum toroidal magnetic field follow from that.
Otherwise we evaluate what are the currents needed for a given flattop duration and toroidal magnetic field.
These currents may or may not exceed the OH and TF current limits.""";
        default=true
    )
    par.j_tolerance = Entry(Real, "", "Tolerance fraction below current limit at which OH and TF operate at"; default=0.4)
    return par
end

"""
    ActorFluxSwing(dd::IMAS.dd, act::ParametersAllActors; kw...)

This actor operate in two ways, depending on `operate_at_j_crit`
* true => Figure out what is the OH and TF current limit, and evaluate flattop duration and maximum toroidal magnetic field follow from that
* false => Evaluate what are the currents needed for a given flattop duration and toroidal magnetic field, which may or may not exceed the OH and TF current limits.

OH flux consumption based on:
* rampup estimate based on Ejima coefficient
* flattop consumption
* vertical field from PF coils

!!! note
    Stores data in `dd.build.flux_swing`, `dd.build.tf`, and `dd.build.oh`

"""
function ActorFluxSwing(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorFluxSwing(kw...)
    actor = ActorFluxSwing(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorFluxSwing(dd::IMAS.dd, par::ParametersActor; kw...)
    logging_actor_init(ActorFluxSwing)
    par = par(kw...)
    return ActorFluxSwing(dd, par, par.operate_at_j_crit, par.j_tolerance)
end

"""
    step(actor::ActorFluxSwing; operate_at_j_crit::Bool=actor.operate_at_j_crit, j_tolerance::Real=actor.j_tolerance, only=:all)

`operate_at_j_crit=true`` makes the OH and TF operate at their current limit (within specified `j_tolerance`).
The flattop duration and toroidal magnetic field follow from that.
Otherwise we evaluate what is the currents needed for a given flattop duration and toroidal magnetic field.
These currents may or may not exceed the OH and TF current limits.

The `only` parameter controls if :tf, :oh, or :all (both) should be calculated
"""
function _step(actor::ActorFluxSwing; operate_at_j_crit::Bool=actor.operate_at_j_crit, j_tolerance::Real=actor.j_tolerance, only=:all)

    bd = actor.dd.build
    target = actor.dd.target
    eq = actor.dd.equilibrium
    eqt = eq.time_slice[]
    cp = actor.dd.core_profiles
    cp1d = cp.profiles_1d[]

    if only ∈ [:all, :oh]
        bd.flux_swing.rampup = rampup_flux_estimates(eqt, cp)
        bd.flux_swing.pf = pf_flux_estimates(eqt)

        if operate_at_j_crit
            oh_maximum_J_B!(bd; j_tolerance)
            bd.flux_swing.flattop = flattop_flux_estimates(bd) # flattop flux based on available current
        else
            bd.flux_swing.flattop = flattop_flux_estimates(target, cp1d) # flattop flux based on target duration
            oh_required_J_B!(bd)
        end

        # flattop duration
        flattop_duration!(bd, eqt, cp1d)
    end

    if only ∈ [:all, :tf]
        tf_required_J_B!(bd, eq)
    end

    return actor
end

"""
    rampup_flux_estimates(eqt::IMAS.equilibrium__time_slice, cp::IMAS.core_profiles)

Estimate OH flux requirement during rampup, where
eqt is supposed to be the equilibrium right at the end of the rampup phase, beginning of flattop
and core_profiles is only used to get core_profiles.global_quantities.ejima
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
    flattop_flux_estimates(target::IMAS.target, cp1d::IMAS.core_profiles__profiles_1d)

Estimate OH flux requirement during flattop (if j_ohmic profile is missing then steady state ohmic profile is assumed)
"""
function flattop_flux_estimates(target::IMAS.target, cp1d::IMAS.core_profiles__profiles_1d)
    return abs(integrate(cp1d.grid.area, cp1d.j_ohmic ./ cp1d.conductivity_parallel)) * target.flattop_duration # V*s
end

"""
    flattop_flux_estimates(bd::IMAS.build; double_swing::Bool=true)

OH flux given its max_b_field and geometry
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

Estimate vertical field from PF coils and its contribution to flux swing, where
`eqt` is supposed to be the equilibrium right at the end of the rampup phase, beginning of flattop
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

#= ========== =#
#  LFS sizing  #
#= ========== =#
mutable struct ActorLFSsizing <: ReactorAbstractActor
    dd::IMAS.dd
    par::ParametersActor
    function ActorLFSsizing(dd::IMAS.dd, par::ParametersActor; kw...)
        logging_actor_init(ActorLFSsizing)
        par = par(kw...)
        return new(dd, par)
    end
end

function ParametersActor(::Type{Val{:ActorLFSsizing}})
    par = ParametersActor(nothing)
    par.do_plot = Entry(Bool, "", "plot"; default=false)
    par.verbose = Entry(Bool, "", "verbose"; default=false)
    return par
end

"""
    ActorLFSsizing(dd::IMAS.dd, act::ParametersAllActors; kw...)

Actor that resizes the Low Field Side of the build.
* Places TF outer leg at radius required to meet the dd.build.tf.ripple requirement
* Other low-field side layers are scaled proportionally

!!! note 
    Manipulates radial build information in `dd.build.layer`
"""
function ActorLFSsizing(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorLFSsizing(kw...)
    actor = ActorLFSsizing(dd, par)
    if par.do_plot
        plot(dd.build)
    end
    step(actor; par.verbose)
    finalize(actor)
    if par.do_plot
        display(plot!(dd.build; cx=false))
    end
    return actor
end

function _step(actor::ActorLFSsizing; verbose::Bool=false)
    dd = actor.dd

    new_TF_radius = IMAS.R_tf_ripple(IMAS.get_build(dd.build, type=_plasma_).end_radius, dd.build.tf.ripple, dd.build.tf.coils_n)

    itf = IMAS.get_build(dd.build, type=_tf_, fs=_lfs_, return_index=true) - 1
    iplasma = IMAS.get_build(dd.build, type=_plasma_, return_index=true) + 1

    # resize layers proportionally
    # start from the vacuum gaps before resizing the material layers
    for vac in [true, false]
        old_TF_radius = IMAS.get_build(dd.build, type=_tf_, fs=_lfs_).start_radius
        delta = new_TF_radius - old_TF_radius
        if verbose
            println("TF radius changed by $delta [m]")
        end
        thicknesses = [dd.build.layer[k].thickness for k in iplasma:itf if !vac || lowercase(dd.build.layer[k].material) == "vacuum"]
        for k in iplasma:itf
            if !vac || lowercase(dd.build.layer[k].material) == "vacuum"
                dd.build.layer[k].thickness *= (1 + delta / sum(thicknesses))
                hfs_thickness = IMAS.get_build(dd.build, identifier=dd.build.layer[k].identifier, fs=_hfs_).thickness
                if dd.build.layer[k].thickness < hfs_thickness
                    dd.build.layer[k].thickness = hfs_thickness
                end
            end
        end
    end

end

#= ========== =#
#  HFS sizing  #
#= ========== =#
mutable struct ActorHFSsizing <: ReactorAbstractActor
    dd::IMAS.dd
    par::ParametersActor
    stresses_actor::ActorStresses
    fluxswing_actor::ActorFluxSwing
end

function ParametersActor(::Type{Val{:ActorHFSsizing}})
    par = ParametersActor(nothing)
    par.j_tolerance = Entry(Float64, "", "Tolerance on the conductor current limits"; default=0.4)
    par.stress_tolerance = Entry(Float64, "", "Tolerance on the structural stresses limits"; default=0.2)
    par.fixed_aspect_ratio = Entry(Bool, "", "Raise an error if aspect_ratio changes more than 10%"; default=true)
    par.unconstrained_flattop_duration = Entry(Bool, "", "Maximize flux_duration without targeting a specific value"; default=true)
    par.do_plot = Entry(Bool, "", "plot"; default=false)
    par.verbose = Entry(Bool, "", "verbose"; default=false)
    return par
end

"""
    ActorHFSsizing(dd::IMAS.dd, act::ParametersAllActors; kw...)

Actor that resizes the High Field Side of the build.
* takes into account the OH maximum allowed superconductor current/Field
* takes into account the stresses on the center stack
    
!!! note 
    Manipulates radial build information in `dd.build.layer`
"""
function ActorHFSsizing(dd::IMAS.dd, act::ParametersAllActors; kw_ActorFluxSwing=Dict(), kw_ActorStresses=Dict(), kw...)
    par = act.ActorHFSsizing(kw...)
    actor = ActorHFSsizing(dd, par, act; kw_ActorFluxSwing, kw_ActorStresses)
    if par.do_plot
        p = plot(dd.build)
    end
    step(actor)
    finalize(actor)
    if par.do_plot
        display(plot!(p, dd.build; cx=false))
    end
    return actor
end

function ActorHFSsizing(dd::IMAS.dd, par::ParametersActor, act::ParametersAllActors; kw_ActorFluxSwing=Dict(), kw_ActorStresses=Dict(), kw...)
    par = act.ActorHFSsizing(kw...)
    fluxswing_actor = ActorFluxSwing(dd, act.ActorFluxSwing; kw_ActorFluxSwing...)
    stresses_actor = ActorStresses(dd, act.ActorStresses; kw_ActorStresses...)
    return ActorHFSsizing(dd, par, stresses_actor, fluxswing_actor)
end

function _step(
    actor::ActorHFSsizing;
    j_tolerance::Real=actor.par.j_tolerance,
    stress_tolerance::Real=actor.par.stress_tolerance,
    fixed_aspect_ratio::Bool=actor.par.fixed_aspect_ratio,
    unconstrained_flattop_duration::Bool=actor.par.unconstrained_flattop_duration,
    verbose::Bool=actor.par.verbose,
    debug_plot=false
)

    function target_value(value, target, tolerance) # relative error with tolerance
        return abs((value .* (1.0 .+ tolerance) .- target) ./ (abs(target) + 1.0))
    end

    function assign_PL_OH_TF_thicknesses(x0, what)
        x0 = map(abs, x0)
        c_extra = 0.0

        if what == :oh
            OH.thickness = x0[1]
            if length(x0) == 2
                fraction_stainless, c_extra = mirror_bound_w_cost(x0[2], 0.5, 1.0 - dd.build.oh.technology.fraction_void - 0.05)
                dd.build.oh.technology.fraction_stainless = fraction_stainless
            end

        elseif what == :tf
            TFhfs.thickness = x0[1]
            if length(x0) == 2
                fraction_stainless, c_extra = mirror_bound_w_cost(x0[2], 0.5, 1.0 - dd.build.oh.technology.fraction_void - 0.05)
                dd.build.tf.technology.fraction_stainless = fraction_stainless
            end

        else
            OH.thickness, TFhfs.thickness = x0
            if length(x0) == 4
                fraction_stainless, c_extra = mirror_bound_w_cost(x0[3], 0.5, 1.0 - dd.build.oh.technology.fraction_void - 0.05)
                dd.build.oh.technology.fraction_stainless = fraction_stainless
                fraction_stainless, c_extra = mirror_bound_w_cost(x0[4], 0.5, 1.0 - dd.build.oh.technology.fraction_void - 0.05)
                dd.build.tf.technology.fraction_stainless = fraction_stainless
            end
        end

        plug.thickness += old_plasma_start_radius - plasma.start_radius
        plug.thickness = max(OH.thickness / 4.0, plug.thickness)

        TFlfs.thickness = TFhfs.thickness
        return c_extra
    end

    function cost(x0, what)
        # assign optimization arguments and evaluate coils currents and stresses
        c_extra = assign_PL_OH_TF_thicknesses(x0, what)
        _step(actor.fluxswing_actor; operate_at_j_crit=unconstrained_flattop_duration, j_tolerance, only=what)
        _step(actor.stresses_actor)

        # OH and plug sizing based on stresses
        c_joh = c_soh = c_spl = 0.0
        if what ∈ [:oh, :all]
            c_joh1 = target_value(dd.build.oh.critical_j, dd.build.oh.max_j, -j_tolerance)
            c_joh2 = target_value(dd.build.oh.max_j, dd.build.oh.critical_j, j_tolerance)
            c_joh = norm([c_joh1, c_joh2])
            c_soh = target_value(maximum(dd.solid_mechanics.center_stack.stress.vonmises.oh), stainless_steel.yield_strength, stress_tolerance)
            if !ismissing(dd.solid_mechanics.center_stack.stress.vonmises, :pl)
                c_spl = target_value(maximum(dd.solid_mechanics.center_stack.stress.vonmises.pl), stainless_steel.yield_strength, stress_tolerance)
            end
        end

        # TF sizing based on stresses
        c_jtf = c_stf = 0.0
        if what ∈ [:tf, :all]
            c_jtf1 = target_value(dd.build.tf.critical_j, dd.build.tf.max_j, -j_tolerance)
            c_jtf2 = target_value(dd.build.tf.max_j, dd.build.tf.critical_j, j_tolerance)
            c_jtf = norm([c_jtf1, c_jtf2])
            c_stf = target_value(maximum(dd.solid_mechanics.center_stack.stress.vonmises.tf), stainless_steel.yield_strength, stress_tolerance)
        end

        if debug_plot
            push!(C_JOH, c_joh)
            push!(C_SOH, c_soh)
            push!(C_JTF, c_jtf)
            push!(C_STF, c_stf)
        end

        # total cost
        return norm(vcat([c_joh, c_jtf], [c_soh, c_stf, c_spl], [c_extra]))
    end

    @assert actor.stresses_actor.dd === actor.fluxswing_actor.dd
    dd = actor.stresses_actor.dd
    target_B0 = maximum(abs.(dd.equilibrium.vacuum_toroidal_field.b0))

    # init
    plug = dd.build.layer[1]
    OH = IMAS.get_build(dd.build, type=_oh_)
    TFhfs = IMAS.get_build(dd.build, type=_tf_, fs=_hfs_)
    TFlfs = IMAS.get_build(dd.build, type=_tf_, fs=_lfs_)
    iplasma = IMAS.get_build(dd.build, type=_plasma_, return_index=true)
    plasma = dd.build.layer[iplasma]

    old_R0 = (TFhfs.end_radius + TFlfs.start_radius) / 2.0
    old_plasma_start_radius = plasma.start_radius
    old_a = plasma.thickness / 2.0
    old_ϵ = old_R0 / old_a

    if debug_plot
        C_JOH = []
        C_SOH = []
        C_JTF = []
        C_STF = []
    end

    # initialize all dd fields
    step(actor.fluxswing_actor; operate_at_j_crit=unconstrained_flattop_duration, j_tolerance)
    step(actor.stresses_actor)

    dd.build.oh.technology.fraction_stainless = 0.5
    dd.build.tf.technology.fraction_stainless = 0.5

    # plug and OH optimization (w/fraction)
    old_thicknesses = [layer.thickness for layer in dd.build.layer]
    res = Optim.optimize(
        x0 -> cost(x0, :oh),
        [OH.thickness, dd.build.oh.technology.fraction_stainless],
        Optim.NelderMead(),
        Optim.Options(time_limit=60);
        autodiff=:forward
    )
    assign_PL_OH_TF_thicknesses(res.minimizer, :oh)
    _step(actor.fluxswing_actor; operate_at_j_crit=unconstrained_flattop_duration, j_tolerance, only=:oh)
    _step(actor.stresses_actor)
    if verbose
        display(res)
    end

    # TF optimization (w/fraction)
    old_thicknesses = [layer.thickness for layer in dd.build.layer]
    res = Optim.optimize(
        x0 -> cost(x0, :tf),
        [TFhfs.thickness, dd.build.tf.technology.fraction_stainless],
        Optim.NelderMead(),
        Optim.Options(time_limit=60);
        autodiff=:forward
    )
    assign_PL_OH_TF_thicknesses(res.minimizer, :tf)
    _step(actor.fluxswing_actor; operate_at_j_crit=unconstrained_flattop_duration, j_tolerance, only=:tf)
    _step(actor.stresses_actor)
    if verbose
        display(res)
    end

    # combined plug+OH+TF optimization
    res = nothing
    if (dd.solid_mechanics.center_stack.bucked == 1 || dd.solid_mechanics.center_stack.noslip == 1 || dd.solid_mechanics.center_stack.plug == 1)
        old_thicknesses = [layer.thickness for layer in dd.build.layer]
        res = Optim.optimize(
            x0 -> cost(x0, :all),
            [OH.thickness, TFhfs.thickness, dd.build.oh.technology.fraction_stainless, dd.build.tf.technology.fraction_stainless],
            Optim.NelderMead(),
            Optim.Options(time_limit=60, iterations=1000);
            autodiff=:forward
        )
        assign_PL_OH_TF_thicknesses(res.minimizer, :all)
        _step(actor.fluxswing_actor; operate_at_j_crit=unconstrained_flattop_duration, j_tolerance)
        _step(actor.stresses_actor)
        if verbose
            display(res)
        end
    end

    R0 = (TFhfs.end_radius + TFlfs.start_radius) / 2.0
    a = plasma.thickness / 2.0
    ϵ = R0 / a

    if debug_plot
        p = plot(yscale=:log10)
        plot!(p, C_JOH ./ (C_JOH .> 0.0), label="Jcrit OH")
        plot!(p, C_SOH ./ (C_SOH .> 0.0), label="Stresses OH")
        plot!(p, C_JTF ./ (C_JTF .> 0.0), label="Jcrit TF")
        plot!(p, C_STF ./ (C_STF .> 0.0), label="Stresses TF")
        display(p)
    end

    if verbose
        R0 = (TFhfs.end_radius + TFlfs.start_radius) / 2.0
        @show target_B0
        @show dd.build.tf.max_b_field * TFhfs.end_radius / R0

        @show dd.build.oh.flattop_duration
        @show dd.target.flattop_duration

        @show dd.build.oh.max_j
        @show dd.build.oh.critical_j

        @show dd.build.tf.max_j
        @show dd.build.tf.critical_j

        @show maximum(dd.solid_mechanics.center_stack.stress.vonmises.oh)
        @show stainless_steel.yield_strength

        @show maximum(dd.solid_mechanics.center_stack.stress.vonmises.tf)
        @show stainless_steel.yield_strength

        @show ϵ
        @show old_ϵ
    end

    function rel_error(value, target) # relative error with tolerance
        return abs((value .- target) ./ target)
    end

    max_B0 = dd.build.tf.max_b_field / TFhfs.end_radius * R0
    @assert target_B0 < max_B0 "TF cannot achieve requested B0 ($target_B0 --> $max_B0)"

    @assert dd.build.oh.max_j < dd.build.oh.critical_j
    @assert dd.build.tf.max_j < dd.build.tf.critical_j
    @assert maximum(dd.solid_mechanics.center_stack.stress.vonmises.oh) < stainless_steel.yield_strength
    @assert maximum(dd.solid_mechanics.center_stack.stress.vonmises.tf) < stainless_steel.yield_strength
    if !unconstrained_flattop_duration
        @assert rel_error(dd.build.oh.flattop_duration, dd.target.flattop_duration) < 0.1 "Relative error on flattop duration is more than 10% ($(dd.build.oh.flattop_duration) --> $(dd.target.flattop_duration))"
    end
    if fixed_aspect_ratio
        @assert rel_error(ϵ, old_ϵ) < 0.1 "ActorHFSsizing: plasma aspect ratio changed more than 10% ($old_ϵ --> $ϵ)"
    end

    return actor
end

#= ============= =#
#  cross-section  #
#= ============= =#
mutable struct ActorCXbuild <: ReactorAbstractActor
    dd::IMAS.dd
    par::ParametersActor
    function ActorCXbuild(dd::IMAS.dd, par::ParametersActor; kw...)
        logging_actor_init(ActorCXbuild)
        par = par(kw...)
        return new(dd, par)
    end
end

function ParametersActor(::Type{Val{:ActorCXbuild}})
    par = ParametersActor(nothing)
    par.rebuild_wall = Entry(Bool, "", "Rebuild wall based on equilibrium"; default=false)
    par.do_plot = Entry(Bool, "", "plot"; default=false)
    return par
end

"""
    ActorCXbuild(dd::IMAS.dd, act::ParametersAllActors; kw...)

Actor that builds the 2D cross section of the build.

!!! note 
    Manipulates data in `dd.build`
"""
function ActorCXbuild(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorCXbuild(kw...)
    actor = ActorCXbuild(dd, par)
    step(actor)
    finalize(actor)
    if par.do_plot
        plot(dd.build)
        display(plot!(dd.build; cx=false))
    end
    return actor
end

function _step(actor::ActorCXbuild; rebuild_wall::Bool=actor.par.rebuild_wall)
    build_cx!(actor.dd; rebuild_wall)
end

"""
    wall_from_eq(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice; divertor_length_fraction::Real=0.2)

Generate first wall outline starting from an equilibrium
"""
function wall_from_eq(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice; divertor_length_fraction::Real=0.2)
    R0 = eqt.global_quantities.magnetic_axis.r
    Z0 = eqt.global_quantities.magnetic_axis.z

    # lcfs
    ψb = IMAS.find_psi_boundary(eqt)
    rlcfs, zlcfs, _ = IMAS.flux_surface(eqt, ψb, true)

    # Set the radial build thickness of the plasma
    plasma = IMAS.get_build(bd, type=_plasma_)
    a = (minimum(rlcfs) - plasma.start_radius)
    plasma.thickness = maximum(rlcfs) - minimum(rlcfs) + 2 * a
    R_hfs_plasma = plasma.start_radius
    R_lfs_plasma = plasma.end_radius

    # main chamber (clip elements that go beyond plasma radial build thickness)
    plasma_poly = xy_polygon(rlcfs, zlcfs)
    wall_poly = LibGEOS.buffer(plasma_poly, a)
    R = [v[1] for v in GeoInterface.coordinates(wall_poly)[1]]
    Z = [v[2] for v in GeoInterface.coordinates(wall_poly)[1]]
    R[R.<R_hfs_plasma] .= R_hfs_plasma
    R[R.>R_lfs_plasma] .= R_lfs_plasma
    Z = (Z .- Z0) .* 1.05 .+ Z0
    wall_poly = xy_polygon(R, Z)

    t = LinRange(0, 2π, 31)

    # divertor lengths
    linear_plasma_size = sqrt((maximum(zlcfs) - minimum(zlcfs)) * (maximum(rlcfs) - minimum(rlcfs)))
    max_divertor_length = linear_plasma_size * divertor_length_fraction

    # private flux regions
    private = IMAS.flux_surface(eqt, ψb, false)
    for (pr, pz) in private
        if sign(pz[1] - Z0) != sign(pz[end] - Z0)
            # open flux surface does not encicle the plasma
            continue
        elseif IMAS.minimum_distance_two_shapes(pr, pz, rlcfs, zlcfs) > linear_plasma_size / 5
            # secondary Xpoint far away
            continue
        end

        # xpoint at location of maximum curvature
        index = argmax(abs.(IMAS.curvature(pr,pz)))
        Rx = pr[index]
        Zx = pz[index]

        max_d = maximum(sqrt.((Rx .- pr) .^ 2.0 .+ (Zx .- pz) .^ 2.0))
        divertor_length = min(max_d * 0.8, max_divertor_length)

        # limit extent of private flux regions
        circle = collect(zip(divertor_length .* cos.(t) .+ Rx, sign(Zx) .* divertor_length .* sin.(t) .+ Zx))
        circle[1] = circle[end]
        slot = [(rr, zz) for (rr, zz) in zip(pr, pz) if PolygonOps.inpolygon((rr, zz), circle) == 1 && rr >= R_hfs_plasma && rr <= R_lfs_plasma]
        pr = [rr for (rr, zz) in slot]
        pz = [zz for (rr, zz) in slot]

        if isempty(pr)
            error("Something is wrong with the geometry and equilibrium")
        end

        # remove private flux region from wall (necessary because of Z expansion)
        wall_poly = LibGEOS.difference(wall_poly, xy_polygon(pr, pz))

        # add the divertor slots
        α = 0.2
        pr = vcat(pr, R0 * α + Rx * (1 - α))
        pz = vcat(pz, Z0 * α + Zx * (1 - α))
        slot = LibGEOS.buffer(xy_polygon(pr, pz), a)
        wall_poly = LibGEOS.union(wall_poly, slot)
    end

    # vertical clip
    wall_poly = LibGEOS.difference(wall_poly, xy_polygon(rectangle_shape(0, R_hfs_plasma, 100)...))
    wall_poly = LibGEOS.difference(wall_poly, xy_polygon(rectangle_shape(R_lfs_plasma, 10 * R_lfs_plasma, 100)...))

    # round corners
    wall_poly = LibGEOS.buffer(wall_poly, -a / 4)
    wall_poly = LibGEOS.buffer(wall_poly, a / 4)

    pr = [v[1] for v in GeoInterface.coordinates(wall_poly)[1]]
    pz = [v[2] for v in GeoInterface.coordinates(wall_poly)[1]]

    pr, pz = IMAS.resample_2d_line(pr, pz; step=0.1)

    return pr, pz
end

function divertor_regions!(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice)
    R0 = eqt.global_quantities.magnetic_axis.r
    Z0 = eqt.global_quantities.magnetic_axis.z

    ipl = IMAS.get_build(bd, type=_plasma_, return_index=true)
    plasma_poly = xy_polygon(bd.layer[ipl])

    wall_poly = xy_polygon(bd.layer[ipl-1])
    for ltype in [_blanket_, _shield_, _wall_,]
        iwl = IMAS.get_build(bd, type=ltype, fs=_hfs_, return_index=true, raise_error_on_missing=false)
        if iwl !== missing
            wall_poly = xy_polygon(bd.layer[iwl])
            break
        end
    end

    divertors = IMAS.IDSvectorElement[]
    for x_point in eqt.boundary.x_point
        Zx = x_point.z
        Rx = x_point.r
        m = (Zx - Z0) / (Rx - R0)
        xx = [0, R0 * 2.0]
        yy = line_through_point(-1.0 ./ m, Rx, Zx, xx)
        pr = vcat(xx, reverse(xx), xx[1])
        pz = vcat(yy, [Zx * 5, Zx * 5], yy[1])

        domain = xy_polygon(pr, pz)
        divertor_poly = LibGEOS.intersection(wall_poly, domain)
        divertor_poly = LibGEOS.difference(divertor_poly, plasma_poly)

        pr = [v[1] for v in GeoInterface.coordinates(divertor_poly)[1]]
        pz = [v[2] for v in GeoInterface.coordinates(divertor_poly)[1]]

        # assign to build structure
        if Zx > Z0
            name = "Upper divertor"
        else
            name = "Lower divertor"
        end
        structure = resize!(bd.structure, "type" => Int(_divertor_), "name" => name)
        structure.material = "Tungsten"
        structure.outline.r = pr
        structure.outline.z = pz
        structure.toroidal_extent = 2pi

        push!(divertors, structure)
    end

    return divertors
end

function blanket_regions!(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice)
    R0 = eqt.global_quantities.magnetic_axis.r

    layers = bd.layer
    iblanket = IMAS.get_build(bd; type=_blanket_, fs=_lfs_, return_index=true, raise_error_on_missing=false)
    if iblanket === missing
        return IMAS.IDSvectorElement[]
    end
    layer = layers[iblanket]
    layer_in = layers[iblanket-1]

    layer_poly = xy_polygon(layer)
    layer_in_poly = xy_polygon(layer_in)
    ring_poly = LibGEOS.difference(layer_poly, layer_in_poly)
    for structure in [structure for structure in bd.structure if structure.type == Int(_divertor_)]
        structure_poly = xy_polygon(structure)
        ring_poly = LibGEOS.difference(ring_poly, structure_poly)
    end

    geometries = LibGEOS.getGeometries(ring_poly)
    blankets = IMAS.IDSvectorElement[]
    for poly in geometries
        coords = GeoInterface.coordinates(poly)
        pr = [v[1] for v in coords[1]]
        pz = [v[2] for v in coords[1]]

        # assign to build structure
        if length(geometries) == 2
            if sum(pr) / length(pr) > R0
                name = "LFS blanket"
            else
                name = "HFS blanket"
            end
        else
            name = "blanket"
        end

        structure = resize!(bd.structure, "type" => Int(_blanket_), "name" => name)
        structure.outline.r = pr
        structure.outline.z = pz
        structure.toroidal_extent = 2pi

        push!(blankets, structure)
    end

    return blankets
end

"""
    build_cx!(dd::IMAS.dd; rebuild_wall::Bool=false)

Translates 1D build to 2D cross-sections starting either wall information
If wall information is missing, then the first wall information is generated starting from equilibrium time_slice
"""
function build_cx!(dd::IMAS.dd; rebuild_wall::Bool=false)
    wall = IMAS.first_wall(dd.wall)
    if wall === missing || rebuild_wall
        pr, pz = wall_from_eq(dd.build, dd.equilibrium.time_slice[])
    else
        pr = wall.r
        pz = wall.z
    end

    build_cx!(dd.build, pr, pz)

    divertor_regions!(dd.build, dd.equilibrium.time_slice[])

    blanket_regions!(dd.build, dd.equilibrium.time_slice[])

    if wall === missing || rebuild_wall
        plasma = IMAS.get_build(dd.build, type=_plasma_)
        resize!(dd.wall.description_2d, 1)
        resize!(dd.wall.description_2d[1].limiter.unit, 1)
        dd.wall.description_2d[1].limiter.unit[1].outline.r = plasma.outline.r
        dd.wall.description_2d[1].limiter.unit[1].outline.z = plasma.outline.z
    end

    return dd.build
end

"""
    build_cx!(bd::IMAS.build, pr::Vector{Float64}, pz::Vector{Float64})

Translates 1D build to 2D cross-sections starting from R and Z coordinates of plasma first wall
"""
function build_cx!(bd::IMAS.build, pr::Vector{Float64}, pz::Vector{Float64})
    ipl = IMAS.get_build(bd, type=_plasma_, return_index=true)
    itf = IMAS.get_build(bd, type=_tf_, fs=_hfs_, return_index=true)

    # _plasma_ outline scaled to match 1D radial build
    start_radius = bd.layer[ipl].start_radius
    end_radius = bd.layer[ipl].end_radius
    pr1 = minimum(pr)
    pr2 = maximum(pr)
    fact = (end_radius - start_radius) / (pr2 - pr1)
    pz .= pz .* fact
    pr .= (pr .- pr1) .* fact .+ start_radius
    bd.layer[ipl].outline.r = pr
    bd.layer[ipl].outline.z = pz

    coils_inside = any([contains(lowercase(l.name), "coils") for l in bd.layer])

    # all layers between plasma and OH
    # k-1 means the layer outside (ie. towards the tf)
    # k   is the current layer
    # k+1 means the layer inside (ie. towards the plasma)
    # 
    # forward pass: from plasma to TF _convex_hull_ and then desired TF shape
    tf_to_plasma = IMAS.get_build(bd, fs=_hfs_, return_only_one=false, return_index=true)
    plasma_to_tf = reverse(tf_to_plasma)
    for k in plasma_to_tf
        original_shape = bd.layer[k].shape
        if k == itf + 1
            # layer that is inside of the TF sets TF shape
            layer_shape = BuildLayerShape(mod(mod(bd.layer[k].shape, 1000), 100))
            optimize_shape(bd, k + 1, k, layer_shape; tight=!coils_inside)
        else
            # everything else is conformal convex hull
            optimize_shape(bd, k + 1, k, _convex_hull_)
            bd.layer[k].shape = original_shape
        end
    end
    # reverse pass: from TF to plasma only with negative offset
    for k in tf_to_plasma[2:end]
        if bd.layer[k+1].shape == Int(_negative_offset_)
            optimize_shape(bd, k, k + 1, _negative_offset_)
        end
    end
    # forward pass: from plasma to TF with desired shapes
    for k in plasma_to_tf[1:end-1]
        if bd.layer[k].shape == Int(_negative_offset_)
            break
        else
            layer_shape = BuildLayerShape(mod(mod(bd.layer[k].shape, 1000), 100))
            optimize_shape(bd, k + 1, k, layer_shape)
        end
    end

    # _in_
    D = minimum(IMAS.get_build(bd, type=_tf_, fs=_hfs_).outline.z)
    U = maximum(IMAS.get_build(bd, type=_tf_, fs=_hfs_).outline.z)
    for k in IMAS.get_build(bd, fs=_in_, return_index=true, return_only_one=false)
        L = bd.layer[k].start_radius
        R = bd.layer[k].end_radius
        bd.layer[k].outline.r, bd.layer[k].outline.z = rectangle_shape(L, R, D, U)
    end

    # _out_
    iout = IMAS.get_build(bd, fs=_out_, return_index=true, return_only_one=false)
    if lowercase(bd.layer[iout[end]].name) == "cryostat"
        olfs = IMAS.get_build(bd, fs=_lfs_, return_index=true, return_only_one=false)[end]
        optimize_shape(bd, olfs, iout[end], _silo_)
        for k in reverse(iout[2:end])
            optimize_shape(bd, k, k - 1, _negative_offset_)
        end
    else
        for k in iout
            L = 0
            R = bd.layer[k].end_radius
            D = minimum(bd.layer[k-1].outline.z) - bd.layer[k].thickness
            U = maximum(bd.layer[k-1].outline.z) + bd.layer[k].thickness
            bd.layer[k].outline.r, bd.layer[k].outline.z = rectangle_shape(L, R, D, U)
        end
    end

    return bd
end

"""
    optimize_shape(bd::IMAS.build, obstr_index::Int, layer_index::Int, shape::BuildLayerShape)

Generates outline of layer in such a way to maintain minimum distance from inner layer
"""
function optimize_shape(bd::IMAS.build, obstr_index::Int, layer_index::Int, shape::BuildLayerShape; tight::Bool=true)
    layer = bd.layer[layer_index]
    obstr = bd.layer[obstr_index]
    # display("Layer $layer_index = $(layer.name)")
    # display("Obstr $obstr_index = $(obstr.name)")
    if layer.fs == Int(_out_)
        l_start = 0
        l_end = layer.end_radius
        o_start = 0
        o_end = obstr.end_radius
    else
        if obstr.fs in [Int(_lhfs_), Int(_out_)]
            o_start = obstr.start_radius
            o_end = obstr.end_radius
        else
            o_start = obstr.start_radius
            o_end = IMAS.get_build(bd, identifier=obstr.identifier, fs=_lfs_).end_radius
        end
        l_start = layer.start_radius
        if layer.type == Int(_plasma_)
            l_end = layer.end_radius
        else
            l_end = IMAS.get_build(bd, identifier=layer.identifier, fs=_lfs_).end_radius
        end
    end
    hfs_thickness = o_start - l_start
    lfs_thickness = l_end - o_end
    oR = obstr.outline.r
    oZ = obstr.outline.z
    if layer.fs == Int(_out_)
        target_minimum_distance = lfs_thickness
    else
        if tight
            target_minimum_distance = min(hfs_thickness, lfs_thickness)
        else
            target_minimum_distance = sqrt(hfs_thickness^2 + lfs_thickness^2) / 2.0
        end
    end
    r_offset = (lfs_thickness .- hfs_thickness) / 2.0

    # update shape
    layer.shape = Int(shape)

    # handle offset, negative offset, negative offset, and convex-hull
    if layer.shape in [Int(_offset_), Int(_negative_offset_), Int(_convex_hull_)]
        poly = LibGEOS.buffer(xy_polygon(oR, oZ), (hfs_thickness + lfs_thickness) / 2.0)
        R = [v[1] .+ r_offset for v in GeoInterface.coordinates(poly)[1]]
        Z = [v[2] for v in GeoInterface.coordinates(poly)[1]]
        if layer.shape == Int(_convex_hull_)
            hull = convex_hull(R, Z; closed_polygon=true)
            R = [r for (r, z) in hull]
            Z = [z for (r, z) in hull]
            # resample disabled because this can lead to outlines of different layers to be crossing
            # R, Z = IMAS.resample_2d_line(R, Z)
        end
        layer.outline.r, layer.outline.z = R, Z

    else # handle shapes
        if layer.shape > 1000
            layer.shape = mod(layer.shape, 1000)
        end
        if layer.shape > 100
            layer.shape = mod(layer.shape, 100)
        end

        if layer.shape == Int(_silo_)
            is_up_down_symmetric = false
        elseif abs(sum(oZ) / sum(abs.(oZ))) < 1E-2
            is_up_down_symmetric = true
        else
            is_up_down_symmetric = false
        end

        is_negative_D = false
        if layer.shape != Int(_silo_)
            _, imaxr = findmax(oR)
            _, iminr = findmin(oR)
            _, imaxz = findmax(oZ)
            _, iminz = findmin(oZ)
            r_at_max_z, max_z = oR[imaxz], oZ[imaxz]
            r_at_min_z, min_z = oR[iminz], oZ[iminz]
            z_at_max_r, max_r = oZ[imaxr], oR[imaxr]
            z_at_min_r, min_r = oZ[iminr], oR[iminr]
            a = 0.5 * (max_r - min_r)
            R = 0.5 * (max_r + min_r)
            δu = (R - r_at_max_z) / a
            δl = (R - r_at_min_z) / a
            if δu + δl < -0.1
                is_negative_D = true
            end
        end

        if is_negative_D
            layer.shape = layer.shape + 1000
        end

        if !is_up_down_symmetric
            layer.shape = layer.shape + 100
        end

        func = shape_function(layer.shape)
        layer.shape_parameters = initialize_shape_parameters(layer.shape, oR, oZ, l_start, l_end, target_minimum_distance)

        layer.outline.r, layer.outline.z = func(l_start, l_end, layer.shape_parameters...)
        layer.shape_parameters = optimize_shape(oR, oZ, target_minimum_distance, func, l_start, l_end, layer.shape_parameters)
        layer.outline.r, layer.outline.z = func(l_start, l_end, layer.shape_parameters...; resample=false)
    end

    IMAS.reorder_flux_surface!(layer.outline.r, layer.outline.z)
    # display(plot!(layer.outline.r, layer.outline.z))
end

function assign_build_layers_materials(dd::IMAS.dd, ini::ParametersAllInits)
    bd = dd.build
    for (k, layer) in enumerate(bd.layer)
        if k == 1 && ini.center_stack.plug
            layer.material = ini.material.wall
        elseif layer.type == Int(_plasma_)
            layer.material = any([layer.type in [Int(_blanket_), Int(_shield_)] for layer in dd.build.layer]) ? "DT_plasma" : "DD_plasma"
        elseif layer.type == Int(_gap_)
            layer.material = "Vacuum"
        elseif layer.type == Int(_oh_)
            layer.material = ini.oh.technology.material
            assign_coil_technology(dd, ini, :oh)
        elseif layer.type == Int(_tf_)
            layer.material = ini.tf.technology.material
            assign_coil_technology(dd, ini, :tf)
        elseif layer.type == Int(_shield_)
            layer.material = ini.material.shield
        elseif layer.type == Int(_blanket_)
            layer.material = ini.material.blanket
        elseif layer.type == Int(_wall_)
            layer.material = ini.material.wall
        elseif layer.type == Int(_vessel_)
            layer.material = "Water, Liquid"
        elseif layer.type == Int(_cryostat_)
            layer.material = ini.material.wall
        end
    end
end

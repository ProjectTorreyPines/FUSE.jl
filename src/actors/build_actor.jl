#= ================ =#
#  flux-swing actor #
#= ================ =#

mutable struct FluxSwingActor <: AbstractActor
    dd::IMAS.dd
end

function FluxSwingActor(dd::IMAS.dd, par::Parameters)
    return FluxSwingActor(dd)
end

# step
function step(flxactor::FluxSwingActor, flattop_duration::Real)
    bd = flxactor.dd.build
    eq = flxactor.dd.equilibrium
    eqt = eq.time_slice[]
    cp = flxactor.dd.core_profiles
    cp1d = cp.profiles_1d[]

    bd.flux_swing_requirements.rampup = rampup_flux_requirements(eqt, cp)
    bd.flux_swing_requirements.flattop = flattop_flux_requirements(cp1d, flattop_duration)
    bd.flux_swing_requirements.pf = pf_flux_requirements(eqt)
    oh_requirements(bd)

    tf_requirements(bd, eq)

    return flxactor
end

"""
    rampup_flux_requirements(eqt::IMAS.equilibrium__time_slice, cp::IMAS.core_profiles)

Estimate OH flux requirement during rampup

NOTES:
* Equations from GASC (Stambaugh FST 2011)
* eqt is supposed to be the equilibrium right at the end of the rampup phase, beginning of flattop
* core_profiles is only used to get core_profiles.global_quantities.ejima
"""
function rampup_flux_requirements(eqt::IMAS.equilibrium__time_slice, cp::IMAS.core_profiles)

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
    flattop_flux_requirements(cp1d::IMAS.core_profiles__profiles_1d, flattop_duration)

Estimate OH flux requirement during flattop 
"""
function flattop_flux_requirements(cp1d::IMAS.core_profiles__profiles_1d, flattop_duration::Real)
    return integrate(cp1d.grid.area, cp1d.j_ohmic ./ cp1d.conductivity_parallel .* flattop_duration) # V*s
end

"""
    pf_flux_requirements(eqt::IMAS.equilibrium__time_slice)

Estimate vertical field from PF coils and its contribution to flux swing

NOTES:
* Equations from GASC (Stambaugh FST 2011)
* eqt is supposed to be the equilibrium right at the end of the rampup phase, beginning of flattop
"""
function pf_flux_requirements(eqt::IMAS.equilibrium__time_slice)
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

"""
    oh_requirements(bd::IMAS.build, double_swing::Bool=true)

Evaluate OH current density and B_field required for rampup and flattop

NOTES:
* Equations from GASC (Stambaugh FST 2011)
* Also relevant: `Engineering design solutions of flux swing with structural requirements for ohmic heating solenoids` Smith, R. A. September 30, 1977
"""
function oh_requirements(bd::IMAS.build, double_swing::Bool = true)
    innerSolenoidRadius = IMAS.get_build(bd, type = _oh_).start_radius
    outerSolenoidRadius = IMAS.get_build(bd, type = _oh_).end_radius
    totalOhFluxReq = bd.flux_swing_requirements.rampup + bd.flux_swing_requirements.flattop + bd.flux_swing_requirements.pf

    # Calculate magnetic field at solenoid bore required to match flux swing request
    RiRoFactor = innerSolenoidRadius / outerSolenoidRadius
    magneticFieldSolenoidBore = 3.0 * totalOhFluxReq / pi / outerSolenoidRadius^2 / (RiRoFactor^2 + RiRoFactor + 1.0) / (double_swing ? 2 : 1)
    currentDensityOH = magneticFieldSolenoidBore / (0.4 * pi * outerSolenoidRadius * (1 - innerSolenoidRadius / outerSolenoidRadius))

    # minimum requirements for OH
    bd.oh.max_b_field = magneticFieldSolenoidBore
    bd.oh.max_j = currentDensityOH * 1E6
    bd.oh.critical_j = coil_Jcrit(bd.oh.max_b_field, bd.oh.technology)
end

"""
    tf_requirements(bd::IMAS.build)

Evaluate TF current density given a B_field
"""
function tf_requirements(bd::IMAS.build, eq::IMAS.equilibrium)
    hfsTF = IMAS.get_build(bd, type = _tf_, fs = _hfs_)
    B0 = abs(maximum(eq.vacuum_toroidal_field.b0))
    R0 = eq.vacuum_toroidal_field.r0

    # current in the TF coils
    current_TF = B0 * R0 * 2pi / constants.μ_0 / bd.tf.coils_n
    TF_cx_area = hfsTF.thickness * bd.tf.wedge_thickness

    bd.tf.max_b_field = B0 * R0 / hfsTF.end_radius
    bd.tf.max_j = current_TF / TF_cx_area
    bd.tf.critical_j = coil_Jcrit(bd.tf.max_b_field, bd.tf.technology)
end

#= ======== =#
#  Stresses  #
#= ======== =#

mutable struct StressesActor <: AbstractActor
    dd::IMAS.dd
end

function step(stressactor::StressesActor; bucked = false, noslip = false, plug = false)
    eq = stressactor.dd.equilibrium
    bd = stressactor.dd.build

    R0 = eq.vacuum_toroidal_field.r0
    B0 = maximum(eq.vacuum_toroidal_field.b0)
    R_tf_in = IMAS.get_build(bd, type = _tf_, fs = _hfs_).start_radius
    R_tf_out = IMAS.get_build(bd, type = _tf_, fs = _hfs_).end_radius
    Bz_oh = bd.oh.max_b_field
    R_oh_in = IMAS.get_build(bd, type = _oh_).start_radius
    R_oh_out = IMAS.get_build(bd, type = _oh_).end_radius
    f_struct_tf = bd.tf.technology.fraction_stainless
    f_struct_oh = bd.oh.technology.fraction_stainless

    smcs = empty!(stressactor.dd.solid_mechanics.center_stack)

    for oh_on in [true, false]
        solve_1D_solid_mechanics!(
            smcs,
            R0,
            B0,
            R_tf_in,
            R_tf_out,
            oh_on ? Bz_oh : 0.0,
            R_oh_in,
            R_oh_out;
            bucked = bucked,
            noslip = noslip,
            plug = plug,
            f_struct_tf = f_struct_tf,
            f_struct_oh = f_struct_oh,
            f_struct_pl = 1.0,
            verbose = false
        )
    end

end

function StressesActor(dd::IMAS.dd, par::Parameters)
    sactor = StressesActor(dd)
    step(sactor; bucked = par.center_stack.bucked, noslip = par.center_stack.noslip, plug = par.center_stack.plug)
    finalize(sactor)
    return sactor
end

@recipe function plot_StressesActor(sactor::StressesActor)
    @series begin
        sactor.dd.solid_mechanics.center_stack.stress
    end
end
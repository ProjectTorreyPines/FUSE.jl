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
function oh_maximum_J_B!(bd::IMAS.build; j_tolerance::Float64)
    OH = IMAS.get_build_layer(bd.layer, type=_oh_)
    innerSolenoidRadius = OH.start_radius
    outerSolenoidRadius = OH.end_radius

    # find maximum superconductor critical_j given self-field
    function max_J_OH(x)
        currentDensityOH = abs(x[1])
        magneticFieldSolenoidBore = currentDensityOH / 1E6 * (0.4 * π * outerSolenoidRadius * (1.0 - innerSolenoidRadius / outerSolenoidRadius))
        critical_j = coil_J_B_crit(magneticFieldSolenoidBore, bd.oh.technology)[1]
        # do not use relative error here. Absolute error tells optimizer to lower currentDensityOH if critical_j==0
        return abs(critical_j - currentDensityOH * (1.0 + j_tolerance))
    end
    res = Optim.optimize(max_J_OH, 0.0, 1E9, Optim.GoldenSection(), rel_tol=1E-3)

    # solenoid maximum current and field
    bd.oh.max_j = abs(res.minimizer[1])
    bd.oh.max_b_field = bd.oh.max_j / 1E6 * (0.4 * π * outerSolenoidRadius * (1.0 - innerSolenoidRadius / outerSolenoidRadius))
    bd.oh.critical_j, bd.oh.critical_b_field = coil_J_B_crit(bd.oh.max_b_field, bd.oh.technology)
    return bd.oh
end

"""
    oh_required_J_B!(bd::IMAS.build; double_swing::Bool=true)

Evaluate OH current density and B_field required to achieve a total flux swing given by: bd.flux_swing.rampup + bd.flux_swing.flattop + bd.flux_swing.pf

NOTES:
* Equations from GASC (Stambaugh FST 2011)
* Also relevant: `Engineering design solutions of flux swing with structural requirements for ohmic heating solenoids` Smith, R. A. September 30, 1977
"""
function oh_required_J_B!(bd::IMAS.build; double_swing::Bool=true)
    OH = IMAS.get_build_layer(bd.layer, type=_oh_)
    innerSolenoidRadius = OH.start_radius
    outerSolenoidRadius = OH.end_radius

    totalOhFluxReq = bd.flux_swing.rampup + bd.flux_swing.flattop + bd.flux_swing.pf

    # Calculate magnetic field at solenoid bore required to match flux swing request
    RiRo_factor = innerSolenoidRadius / outerSolenoidRadius
    magneticFieldSolenoidBore = 3.0 * totalOhFluxReq / π / outerSolenoidRadius^2 / (RiRo_factor^2 + RiRo_factor + 1.0) / (double_swing ? 2.0 : 1.0)
    currentDensityOH = magneticFieldSolenoidBore / (0.4 * π * outerSolenoidRadius * (1.0 - innerSolenoidRadius / outerSolenoidRadius))

    # Minimum requirements for OH
    bd.oh.max_b_field = magneticFieldSolenoidBore
    bd.oh.max_j = currentDensityOH * 1E6
    bd.oh.critical_j, bd.oh.critical_b_field = coil_J_B_crit(bd.oh.max_b_field, bd.oh.technology)
    return bd.oh
end

"""
    flattop_duration!(bd::IMAS.build, cp1d::IMAS.core_profiles__profiles_1d; double_swing::Bool=true)

Estimate OH flux requirement during flattop

NOTE:
* if j_ohmic profile is missing then steady state ohmic profile is assumed
"""
function flattop_duration!(bd::IMAS.build, cp1d::IMAS.core_profiles__profiles_1d; double_swing::Bool=true)
    OH = IMAS.get_build_layer(bd.layer, type=_oh_)
    innerSolenoidRadius = OH.start_radius
    outerSolenoidRadius = OH.end_radius

    # estimate oh flattop flux and duration
    RiRo_factor = innerSolenoidRadius / outerSolenoidRadius
    totalOhFlux = bd.oh.max_b_field * (π * outerSolenoidRadius^2 * (RiRo_factor^2 + RiRo_factor + 1.0) * (double_swing ? 2 : 1)) / 3.0
    bd.flux_swing.flattop = totalOhFlux - bd.flux_swing.rampup - bd.flux_swing.pf
    bd.oh.flattop_duration = bd.flux_swing.flattop / abs(integrate(cp1d.grid.area, cp1d.j_ohmic ./ cp1d.conductivity_parallel))
    return bd.oh
end
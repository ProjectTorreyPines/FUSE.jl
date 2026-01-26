#= ========= =#
#  TF magnet  #
#= ========= =#
"""
    tf_maximum_J_B!(bd::IMAS.build; coil_j_margin::Float64)

Evaluate maxium TF current density and magnetic field for given geometry and technology.

NOTE: This function is here but it's not really used for machine desing, since generally the maximum B0 is a input design parameter.
"""
function tf_maximum_J_B!(bd::IMAS.build; coil_j_margin::Float64)
    hfsTF = IMAS.get_build_layer(bd.layer; type=_tf_, fs=_hfs_)
    conductor_start_radius = hfsTF.start_radius + bd.tf.nose_hfs_fraction * hfsTF.thickness
    TF_cx_area = π * (hfsTF.end_radius^2 - conductor_start_radius^2) / bd.tf.coils_n
    mat_tf = Material(bd.tf.technology)

    # find maximum superconductor critical_j given self-field
    function max_J_TF(x)
        currentDensityTF = abs(x[1])
        current_TF = currentDensityTF * TF_cx_area
        max_b_field = current_TF / hfsTF.end_radius / 2π * IMAS.mks.μ_0 * bd.tf.coils_n
        critical_j = mat_tf.critical_current_density(; Bext=max_b_field)
        # do not use relative error here. Absolute error tells optimizer to lower currentDensityTF if critical_j==0
        return abs(critical_j - currentDensityTF * (1.0 + coil_j_margin))
    end
    res = Optim.optimize(max_J_TF, 0.0, 1E9, Optim.Brent(); rel_tol=1E-3)

    # tf maximum current and field
    bd.tf.max_j = abs(res.minimizer[1])
    current_TF = bd.tf.max_j * TF_cx_area
    bd.tf.max_b_field = current_TF / hfsTF.end_radius / 2π * IMAS.mks.μ_0 * bd.tf.coils_n
    bd.tf.critical_j = mat_tf.critical_current_density(; Bext=bd.tf.max_b_field)
    bd.tf.critical_b_field = mat_tf.critical_magnetic_field(; Bext=bd.tf.max_b_field)
    return bd.tf
end

"""
    tf_required_J_B!(bd::IMAS.build, eq::IMAS.equilibrium)

Evaluate TF current density needed to obtain the maximum required B0 at R0
"""
function tf_required_J_B!(bd::IMAS.build, eq::IMAS.equilibrium)
    hfsTF = IMAS.get_build_layer(bd.layer; type=_tf_, fs=_hfs_)
    plasma = IMAS.get_build_layer(bd.layer; type=_plasma_)
    R0 = (plasma.end_radius + plasma.start_radius) / 2.0
    B0 = maximum(abs, eq.vacuum_toroidal_field.b0)
    mat_tf = Material(bd.tf.technology)

    # current in the TF coils
    current_TF = B0 * R0 * 2π / IMAS.mks.μ_0 / bd.tf.coils_n
    conductor_start_radius = hfsTF.start_radius + bd.tf.nose_hfs_fraction * hfsTF.thickness
    TF_cx_area = π * (hfsTF.end_radius^2 - conductor_start_radius^2) / bd.tf.coils_n

    bd.tf.max_b_field = B0 * R0 / hfsTF.end_radius
    bd.tf.max_j = current_TF / TF_cx_area
    bd.tf.critical_j = mat_tf.critical_current_density(; Bext=bd.tf.max_b_field)
    bd.tf.critical_b_field = mat_tf.critical_magnetic_field(; Bext=bd.tf.max_b_field)
    return bd.tf
end

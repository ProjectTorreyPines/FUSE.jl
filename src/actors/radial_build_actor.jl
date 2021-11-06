
"""
    oh_actor(dd::IMAS.dd, time::Real=0.0; ejima, flux_multiplier=1.0, lswing=2)

Evaluate flux consumption and corresponding max magnetic field at the center of the OH coil solenoid
"""
function oh_actor(dd::IMAS.dd, time::Real=0.0; ejima, flux_multiplier=1.0, lswing=2)

    # from IMAS dd to local variables
    time_index = get_time_index(dd.equilibrium.time_slice, time)
    majorRadius = dd.equilibrium.time_slice[time_index].boundary.geometric_axis.r
    minorRadius = dd.equilibrium.time_slice[time_index].boundary.minor_radius
    elongation = dd.equilibrium.time_slice[time_index].boundary.elongation
    plasmaCurrent = dd.equilibrium.time_slice[time_index].global_quantities.ip / 1E6 # in [MA]
    betaP = dd.equilibrium.time_slice[time_index].global_quantities.beta_pol
    li = dd.equilibrium.time_slice[time_index].global_quantities.li_3 # what li ?
    innerSolenoidRadius, outerSolenoidRadius = IMAS.radial_build_start_end_radii(dd.radial_build, 1)
    if ejima === nothing
        ejima = dd.core_profiles.global_quantities.ejima[time_index]
    end

    # ============================= #
    # evaluate plasma inductance
    plasmaInductanceInternal = 0.4 * 0.5 * pi * majorRadius * li
    plasmaInductanceExternal = 0.4 * pi * majorRadius * (log(8.0 * majorRadius / minorRadius / sqrt(elongation)) - 2.0)
    plasmaInductanceTotal = plasmaInductanceInternal + plasmaInductanceExternal

    # evaluate vertical field and its contribution to flux swing
    verticalFieldAtCenter = 0.1 * plasmaCurrent / majorRadius * (log(8.0 * majorRadius / (minorRadius * sqrt(elongation))) - 1.5 + betaP + 0.5 * li)
    fluxFromVerticalField = 0.8 * verticalFieldAtCenter * pi * (majorRadius^2 - (majorRadius - minorRadius)^2)

    # flux swing
    rampUpFlux = (ejima * 0.4 * pi * majorRadius + plasmaInductanceTotal) * plasmaCurrent

    # required flux swing from OH
    totalRampUpFluxReq = rampUpFlux * flux_multiplier - fluxFromVerticalField

    # Calculate magnetic field at solenoid bore required to match flux swing request (Number of swings in OH coil [1=single, 2=double swing])
    RiRoFactor = innerSolenoidRadius / outerSolenoidRadius
    magneticFieldSolenoidBore = 3.0 * totalRampUpFluxReq / pi / outerSolenoidRadius^2 / (RiRoFactor^2 + RiRoFactor + 1.0) / lswing
    # ============================= #

    # assign max B OH to radial_build IDS
    dd.radial_build.oh_b_field_max = magneticFieldSolenoidBore
    return dd
end


function stress_calculations(dd::IMAS.dd)
    error("not completed yet")
    B0_TF = dd.radial_build.tf_b_field_max
    R0_TF = sum(IMAS.radial_build_start_end_radii(dd.radial_build, -1)) / 2.0
    (Rtf1, Rtf2) = IMAS.radial_build_start_end_radii(dd.radial_build, 2)
    B0_OH = dd.radial_build.oh_b_field_max
    (R_sol1, R_sol2) = IMAS.radial_build_start_end_radii(dd.radial_build, 1)
    s_ax_ave = something
    f_t_ss_tot_in = something
    f_oh_cu_in = something
    f_oh_sa_sh_in = something
    ibuck = something
    stress_calculations(B0_TF, R0_TF, Rtf1, Rtf2, B0_OH, R_sol1, R_sol2, s_ax_ave, f_t_ss_tot_in, f_oh_cu_in, f_oh_sa_sh_in, ibuck)
end


function stress_calculations(
    B0_TF, # magnetic field on axis
    R0_TF, # major radius
    Rtf1,  # inner radius TF
    Rtf2,  # outer radius TF
    B0_OH, # magnetic field solenoid bore
    R_sol1,  # inner solenoid radius
    R_sol2,  # outer solenoid radius
    s_ax_ave,   # average stress axial TF
    f_t_ss_tot_in, # fraction copper + fraction stainless TF
    f_oh_cu_in, # fraction copper + fraction stainless OH
    f_oh_sa_sh_in, # 0.37337
    ibuck) # has plug

    plug_switch = 1
    if ibuck > 1
        plug_switch = ibuck - 1
    end

    robo_tf = B0_TF * R0_TF
    mu0 = 4 * pi * 0.0000001
    r_2 = 0.5 * (Rtf1 + R_sol2)
    r_3 = Rtf2
    em_tf = 193103448275.862
    g_tf = 0.33
    s_t_hoop_ave = -2 / 3 * robo_tf^2 * (2 * r_2 + r_3) / (mu0 * (r_3 - r_2) * (r_3 + r_2)^2)
    f_t_ax_hoop = -s_ax_ave / s_t_hoop_ave
    area_t_ax = pi * (Rtf2^2 - Rtf1^2)
    f_t_ax = s_ax_ave * area_t_ax
    sw_sip1_noslp2 = 1

    b_cs = [0.,B0_OH]

    Rcs_i = R_sol1
    Rcs_o = r_2

    s_c_hoop_ave = b_cs^2 / 6 / mu0 * (Rcs_o + 2 * Rcs_i) / (Rcs_o - Rcs_i)
    f_c_ax_hoop = f_oh_sa_sh_in
    s_c_ax_ave = -f_c_ax_hoop * s_c_hoop_ave
    area_c_ax = pi * (Rcs_o^2 - Rcs_i^2)
    f_c_ax = s_c_ax_ave * area_c_ax

    em_tf = em_tf
    g_tf = g_tf
    s_p_ax_ave = 0
    area_p_ax = pi * (Rcs_i^2)
    f_p_ax = s_p_ax_ave * area_p_ax

    if sw_sip1_noslp2 <= 1
        area_t_ax_use = area_t_ax
    else
        if plug_switch <= 1
            area_t_ax_use = area_t_ax + area_c_ax
        else
            area_t_ax_use = area_t_ax + area_c_ax + area_p_ax
        end
    end
    if sw_sip1_noslp2 <= 1
        f_t_ax_use = f_t_ax
    else
        if plug_switch <= 1
            f_t_ax_use = f_t_ax + f_c_ax
        else
            f_t_ax_use = f_t_ax + f_c_ax + f_p_ax
        end
    end

    s_t_ax_use_nov = f_t_ax_use / area_t_ax_use
    f_t_ss_tot = f_t_ss_tot_in
    s_t_ax_void = s_t_ax_use_nov / f_t_ss_tot
    sw_cs_use = 0

    if sw_sip1_noslp2 <= 1
        area_c_ax_use = area_c_ax
    else
        if plug_switch <= 1
            area_c_ax_use = area_t_ax + area_c_ax
        else
            area_c_ax_use = area_t_ax + area_c_ax + area_p_ax
        end
    end

    if sw_sip1_noslp2 <= 1
        f_c_ax_use = f_c_ax
    else
        if plug_switch <= 1
            f_c_ax_use = f_t_ax + f_c_ax
        else
            f_c_ax_use = f_t_ax + f_c_ax + f_p_ax
        end
    end
    s_c_ax_use_nov = f_c_ax_use / area_c_ax_use

    frac_c_ss_tot = f_oh_cu_in
    s_c_ax_void = s_c_ax_use_nov / frac_c_ss_tot

    if sw_sip1_noslp2 <= 1
        area_p_ax_use = area_p_ax
    else
        if plug_switch <= 1
            area_p_ax_use = area_t_ax + area_c_ax
        else
            area_p_ax_use = area_t_ax + area_c_ax + area_p_ax
        end
    end

    if sw_sip1_noslp2 <= 1
        f_p_ax_use = f_p_ax
    else
        if plug_switch <= 1
            f_p_ax_use = f_t_ax + f_c_ax
        else
            f_p_ax_use = f_t_ax + f_c_ax + f_p_ax
        end
    end

    s_p_ax_use_nov = f_p_ax_use / area_p_ax_use
    frac_p_ss_tot = 1
    s_p_ax_void = s_p_ax_use_nov / frac_p_ss_tot
    C_T = 2 * (1 - g_tf^2) * (R0_TF * B0_TF)^2 / (mu0 * em_tf * (r_3^2 - r_2^2)^2)
    C_C = -(1 - g_tf^2) * b_cs^2 / mu0 / em_tf / (Rcs_o - Rcs_i)^2
    C_P = C_C
    Ebar_tf = em_tf / (1 - g_tf^2)
    Ebar_cs = em_tf / (1 - g_tf^2)
    Ebar_pl = em_tf / (1 - g_tf^2)
    Ebar_cp = Ebar_cs / Ebar_pl / (1 + g_tf)
    Cts3 = Ebar_tf * C_T * ( (3 + g_tf) / 8 * r_3^2 - r_2^2 / 2 * ( (1 + g_tf) * log(r_3) + (1 - g_tf) / 2 ) )
    Ats3 = Ebar_tf * (1 + g_tf)
    Bts3 = -Ebar_tf * (1 - g_tf) / r_3^2
    Cbar_ts3 = -Cts3 / Bts3
    Abar_ts3 = -Ats3 / Bts3
    Ctu2 = C_T * ( r_2^2 / 8 - r_2^2 * ( log(r_2) / 2 - 1.0 / 4.0) )
    Atu2 = 1
    Btu2 = 1 / r_2^2
    Cts2 = Ebar_tf * C_T * ( (3 + g_tf) / 8 * r_2^2 - r_2^2 / 2 * ( (1 + g_tf) * log(r_2) + (1 - g_tf) / 2 ) )
    Ats2 = Ebar_tf * (1 + g_tf)
    Bts2 = -Ebar_tf * (1 - g_tf) / r_2^2
    Atu = Atu2 + Btu2 * Abar_ts3
    Ats = Ats2 + Bts2 * Abar_ts3
    Ccs1 = Ebar_cs * C_C * (Rcs_o * Rcs_i / 3 * (2 + g_tf) - Rcs_i^2 / 8 * (3 + g_tf))
    Acs1 = Ebar_cs * (1 + g_tf)
    Bcs1 = -Ebar_cs * (1 - g_tf) / Rcs_i^2
    Cbar_cs1 = -Ccs1 / Bcs1
    Abar_cs1 = -Acs1 / Bcs1
    CC1 = C_C * ( Ebar_cp * ( (2 + g_tf) / 3 * Rcs_i * Rcs_o - (3 + g_tf) / 8 * Rcs_i^2 ) - (Rcs_i * Rcs_o / 3 - Rcs_i^2 / 8) )
    Ac1 = Ebar_cp * (1 + g_tf) - 1
    Bc1 = -Ebar_cp * (1 - g_tf) / Rcs_i^2 - 1 / Rcs_i^2
    Cbar_c1 = -CC1 / Bc1
    Abar_c1 = -Ac1 / Bc1
    if plug_switch# == 2
        Cbar_c1_use = Cbar_c1
        Abar_c1_use = Abar_c1
    else
        Cbar_c1_use = Cbar_cs1
        Abar_c1_use = Abar_cs1
    end

    Ccu2 = C_C * ( Rcs_o^2 / 3 - Rcs_o^2 / 8)
    Acu2 = 1
    Bcu2 = 1 / Rcs_o^2
    Ccs2 = Ebar_cs * C_C * (Rcs_o * Rcs_o / 3 * (2 + g_tf) - Rcs_o^2 / 8 * (3 + g_tf))
    Acs2 = Ebar_cs * (1 + g_tf)
    Bcs2 = -Ebar_cs * (1 - g_tf) / Rcs_o^2
    Cu = Ctu2 + Btu2 * Cbar_ts3 - (Ccu2 + Bcu2 * Cbar_c1_use)
    Cs = Cts2 + Bts2 * Cbar_ts3 - (Ccs2 + Bcs2 * Cbar_c1_use)
    Acu = Acu2 + Bcu2 * Abar_c1_use
    Acs = Acs2 + Bcs2 * Abar_c1_use
    A_T = (Acu * Cs - Acs * Cu) / (Acs * Atu - Acu * Ats)
    B_T = Cbar_ts3 + Abar_ts3 * A_T
    A_C = (Atu * Cs - Ats * Cu) / (Acs * Atu - Acu * Ats)
    B_C = Cbar_c1_use + Abar_c1_use * A_C

    if plug_switch# == 2:
        A_P = C_C * (Rcs_o * Rcs_i / 3 - Rcs_i^2 / 8) + A_C + B_C / Rcs_i^2
    else
        A_P = 0.
    end
    B_P = 0
    R_min_t = r_2
    u_r_rmin_t = C_T * (R_min_t^2 / 8 - r_2^2 / 2 * (log(R_min_t) - 0.5) )  + A_T + B_T / R_min_t^2
    du_dr_rmin_t = C_T * (3 * R_min_t^2 / 8 - r_2^2 / 2 * (log(R_min_t) + 0.5) )  + A_T - B_T / R_min_t^2
    sr_rmin_t = em_tf / (1 - g_tf^2) * (C_T * (    (3 + g_tf) / 8 * R_min_t^2   - r_2^2 / 2 * ( log(R_min_t) * (1 + g_tf) + (1 - g_tf) / 2 ) )  + A_T * (1 + g_tf) - B_T * (1 - g_tf) / R_min_t^2)
    sh_rmin_t = em_tf / (1 - g_tf^2) * (C_T * (  (1 + 3 * g_tf) / 8 * R_min_t^2 - r_2^2 / 2 * ( log(R_min_t) * (1 + g_tf) - (1 - g_tf) / 2 ) )   + A_T * (1 + g_tf) + B_T * (1 - g_tf) / R_min_t^2)
    svm_t = np.sqrt(((sh_rmin_t - s_t_ax_use_nov)^2 + (s_t_ax_use_nov - sr_rmin_t)^2 + (sr_rmin_t - sh_rmin_t)^2) / 2)
    svm_vd_t = svm_t / f_t_ss_tot
    svm_vd_mp_t = svm_vd_t * 0.000001
    svm_vd_ksi_t = svm_vd_t * 0.000000145
    R_min_c = Rcs_i
    u_r_min_c = C_C * ( Rcs_o * R_min_c / 3 - R_min_c^2 / 8 )  + A_C + B_C / R_min_c^2
    du_dr_rmin_c = C_C * ( 2 * Rcs_o * R_min_c / 3 - 3 * R_min_c^2 / 8 )  + A_C - B_C / R_min_c^2
    sr_rmin_c = em_tf / (1 - g_tf^2) * ( C_C * (   Rcs_o * R_min_c / 3 * (2 + g_tf) - (3 + g_tf) / 8 * R_min_c^2   ) + A_C * (1 + g_tf) - B_C * (1 - g_tf) / R_min_c^2)
    sh_rmin_c = em_tf / (1 - g_tf^2) * ( C_C * (    Rcs_o * R_min_c / 3 * (1 + 2 * g_tf) - (1 + 3 * g_tf) / 8 * R_min_c^2  ) + A_C * (1 + g_tf) + B_C * (1 - g_tf) / R_min_c^2)
    svm_c = np.sqrt(((sh_rmin_c - s_c_ax_use_nov)^2 + (s_c_ax_use_nov - sr_rmin_c)^2 + (sr_rmin_c - sh_rmin_c)^2) / 2)
    svm_vd_c = svm_c / frac_c_ss_tot
    svm_vd_mp_c = svm_vd_c * 0.000001
    svm_vd_ksi_c = svm_vd_c * 0.000000145
#    print (svm_vd_ksi_c)
    R_min_p = 0
    u_r_rmin_p = A_P
    du_dr_rmin_p = A_P
    sr_rmin_p = em_tf * (1 + g_tf) / (1 - g_tf^2) * A_P
    sh_rmin_p = sr_rmin_p
    svm_p = np.sqrt(((sh_rmin_p - s_p_ax_use_nov)^2 + (s_p_ax_use_nov - sr_rmin_p)^2 + (sr_rmin_p - sh_rmin_p)^2) / 2)
    svm_vd_p = svm_p / frac_p_ss_tot
    svm_vd_mp_p = svm_vd_p * 0.000001
    svm_vd_ksi_p = svm_vd_p * 0.000000145

    vals = Dict()
    vals["TF Hoop Stress"] = maximum(sh_rmin_t)
    vals["TF Fraction SS"] = maximum(f_t_ss_tot)
    vals["TF Von Mises Stress"] = maximum(svm_vd_t)
    vals["OH Von Mises Stress"] = maximum(svm_vd_c)
    vals["Plug Von Mises Stress"] = maximum(svm_vd_p)
    vals["OH Buck Switch"] = 0

    return (vals)
end
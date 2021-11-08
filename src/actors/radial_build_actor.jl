using LibGEOS
using Interpolations
using Contour

#= ==== =#
#  init  #
#= ==== =#

"""
    init(radial_build::IMAS.radial_build; Bmax_OH=nothing, Bmax_TF=nothing, layers...)

Initialize radial_build IDS based on center stack layers (thicknesses) and maximum fields

NOTE: index and material follows from standard naming of layers
*   0 ...gap... : vacuum
*   1 OH: ohmic coil
*  -2 inner_TF: toroidal field coil
*  -3 inner_shield...: neutron shield
*  -4 inner_blanket...: neutron blanket
*  -5 inner_wall....: 
*  -1 ...vessel...: 
*   5 outer_wall....: 
*   4 outer_blanket...: neutron blanket
*   3 outer_shield...: neutron shield
*   2 outer_TF: toroidal field coil
"""
function init(radial_build::IMAS.radial_build; Bmax_OH=nothing, Bmax_TF=nothing, layers...)
    if Bmax_OH !== nothing
        radial_build.oh_b_field_max = Bmax_OH
    end
    if Bmax_TF !== nothing
        radial_build.tf_b_field_max = Bmax_TF
    end

    # assign layers
    resize!(radial_build.center_stack, length(layers))
    for (klayer, (layer_name, layer_thickness)) in enumerate(layers)
        radial_build.center_stack[klayer].thickness = layer_thickness
        radial_build.center_stack[klayer].name = replace(String(layer_name), "_" => " ")
        if occursin("gap", lowercase(radial_build.center_stack[klayer].name))
            radial_build.center_stack[klayer].index = 0
            radial_build.center_stack[klayer].material = "vacuum"
        elseif uppercase(radial_build.center_stack[klayer].name) == "OH"
            radial_build.center_stack[klayer].index = 1
        elseif occursin("TF", uppercase(radial_build.center_stack[klayer].name))
            radial_build.center_stack[klayer].index = 2
        elseif occursin("shield", lowercase(radial_build.center_stack[klayer].name))
            radial_build.center_stack[klayer].index = 3
        elseif occursin("blanket", lowercase(radial_build.center_stack[klayer].name))
            radial_build.center_stack[klayer].index = 4
        elseif occursin("wall", lowercase(radial_build.center_stack[klayer].name))
            radial_build.center_stack[klayer].index = 5
        end
        if occursin("inner", lowercase(radial_build.center_stack[klayer].name))
            radial_build.center_stack[klayer].index *= -1
        end
        if occursin("vessel", lowercase(radial_build.center_stack[klayer].name))
            radial_build.center_stack[klayer].index = -1
            radial_build.center_stack[klayer].material = "vacuum"
        end
    end

    return radial_build
end

"""
    init(radial_build::IMAS.radial_build, eqt::IMAS.equilibrium__time_slice; is_nuclear_facility=true)

Simple initialization of radial_build IDS based on equilibrium time_slice
"""
function init(radial_build::IMAS.radial_build, eqt::IMAS.equilibrium__time_slice; is_nuclear_facility=true)
    Bmax_OH = nothing
    Bmax_TF = eqt.profiles_1d.f[end] / eqt.boundary.geometric_axis.r

    rmin = eqt.boundary.geometric_axis.r - eqt.boundary.minor_radius
    rmax = eqt.boundary.geometric_axis.r + eqt.boundary.minor_radius

    if is_nuclear_facility
        n_inner_layers = 6
        gap = (rmax - rmin) / 10.0
        rmin -= gap
        rmax += gap
        dr = rmin / n_inner_layers
        init(radial_build,
            Bmax_OH=Bmax_OH,
            Bmax_TF=Bmax_TF,
            gap_TF=dr * 2.0,
            OH=dr,
            inner_TF=dr,
            inner_shield=dr / 2.0,
            inner_blanket=dr,
            inner_wall=dr / 2.0,
            vacuum_vessel=rmax - rmin,
            outer_wall=dr / 2.0,
            outer_blanket=dr,
            outer_shield=dr / 2.0,
            outer_TF=dr,
            gap_cryostat=3 * dr)

    else
        n_inner_layers = 4.5
        gap = (rmax - rmin) / 10.0
        rmin -= gap
        rmax += gap
        dr = rmin / n_inner_layers
        init(radial_build,
            Bmax_OH=Bmax_OH,
            Bmax_TF=Bmax_TF,
            gap_TF=dr * 2.0,
            OH=dr,
            inner_TF=dr,
            inner_wall=dr / 2.0,
            vacuum_vessel=rmax - rmin,
            outer_wall=dr / 2.0,
            outer_TF=dr,
            gap_cryostat=2 * dr)
    end

    radial_build_cx(radial_build, eqt)

    return radial_build
end

function xy_polygon(x,y)
    if x[1]≈x[end]
        x[end]=x[1]
        y[end]=y[1]
    elseif x[1]!=x[end]
        push!(x,x[1])
        push!(y,y[1])
    end
    coords = [collect(map(collect,zip(x,y)))]
    return LibGEOS.Polygon(coords)
end

function miller(R0,epsilon,kappa,delta,n)
    θ = range(0,2pi,length=n)
    δ₀ = asin(delta)
    x = R0*(1 .+ epsilon .* cos.(θ .+ δ₀*sin.(θ)))
    y = R0*(epsilon*kappa*sin.(θ))
    return (x,y)
end

function wall_miller_conformal(rb,layer_id,elongation,triangularity)
    layer_id=-abs(layer_id)
    Rstart_outer=IMAS.get_radial_build(rb, -layer_id).start_radius
    Rend_outer=IMAS.get_radial_build(rb, -layer_id).end_radius
    if layer_id==-1
        outer_line=miller((Rend_outer+Rstart_outer)/2.0,(Rend_outer-Rstart_outer)/(Rend_outer+Rstart_outer),elongation,triangularity,100)        
        return outer_line,outer_line
    else
        Rstart_inner=IMAS.get_radial_build(rb, layer_id).start_radius
        Rend_inner=IMAS.get_radial_build(rb, layer_id).end_radius
        inner_line=miller((Rstart_outer+Rend_inner)/2.0,(Rstart_outer-Rend_inner)/(Rstart_outer+Rend_inner),elongation,triangularity,100)
        outer_line=miller((Rend_outer+Rstart_inner)/2.0,(Rend_outer-Rstart_inner)/(Rend_outer+Rstart_inner),elongation,triangularity,100)
        return inner_line,outer_line
    end
end

function wall_plug(rb::IMAS.radial_build)
    L=0
    R=IMAS.get_radial_build(rb, 1).start_radius
    U=maximum(IMAS.get_radial_build(rb, -2).outline.z)
    D=minimum(IMAS.get_radial_build(rb, -2).outline.z)
    return [L,R,R,L,L],[D,D,U,U,D]
end

function wall_oh(rb::IMAS.radial_build)
    L=IMAS.get_radial_build(rb, 1).start_radius
    R=IMAS.get_radial_build(rb, 1).end_radius
    U=maximum(IMAS.get_radial_build(rb, -2).outline.z)
    D=minimum(IMAS.get_radial_build(rb, -2).outline.z)
    return [L,R,R,L,L],[D,D,U,U,D]
end

function wall_cryostat(rb::IMAS.radial_build)
    L=0
    R=rb.center_stack[end].end_radius
    U=maximum(IMAS.get_radial_build(rb, -2).outline.z)+rb.center_stack[end].thickness
    D=minimum(IMAS.get_radial_build(rb, -2).outline.z)-rb.center_stack[end].thickness
    return [L,R,R,L,L],[D,D,U,U,D]
end


function radial_build_cx(rb::IMAS.radial_build,eqt,δψ=0.05)
    # we make the outer wall to be conformal to miller
    inner_wall_line, outer_wall_line=wall_miller_conformal(rb,5,eqt.boundary.elongation,eqt.boundary.triangularity) # wall
    wall_poly = xy_polygon(outer_wall_line...)
    if false
        vessel_poly = LibGEOS.buffer(wall_poly,-IMAS.get_radial_build(rb, 5).thickness)
    else
        r = range(eqt.profiles_2d[1].grid.dim1[1], eqt.profiles_2d[1].grid.dim1[end], length=length(eqt.profiles_2d[1].grid.dim1))
        z = range(eqt.profiles_2d[1].grid.dim2[1], eqt.profiles_2d[1].grid.dim2[end], length=length(eqt.profiles_2d[1].grid.dim2))
        PSI_interpolant = Interpolations.CubicSplineInterpolation((r, z), eqt.profiles_2d[1].psi)

        # Inner/outer radii of the vacuum vessel
        R_inner_vessel=IMAS.get_radial_build(rb, -1).start_radius
        R_outer_vessel=IMAS.get_radial_build(rb, -1).end_radius
        psi_vessel_trace=(eqt.global_quantities.psi_boundary-eqt.global_quantities.psi_axis)*(1+δψ)+eqt.global_quantities.psi_axis

        # Trace contours of psi on the inner radii and use those contours as the shape of the vacuum vessel.
        cl = Contour.contour(r,z,eqt.profiles_2d[1].psi,psi_vessel_trace)
        distances=[]
        for line in Contour.lines(cl)
            pr, pz = Contour.coordinates(line)
            push!(distances,minimum(sqrt.((pr.-R_inner_vessel).^2 + (pz.-0.0).^2)))
        end
        inner_vessel_line = Contour.coordinates(Contour.lines(cl)[argmin(distances)])
        outer_vessel_line = Contour.coordinates(Contour.lines(cl)[argmax(distances)])
        inner_vessel_line[1].+=(R_inner_vessel-inner_vessel_line[1][argmin(abs.(inner_vessel_line[2]))])
        outer_vessel_line[1].+=(R_outer_vessel-outer_vessel_line[1][argmin(abs.(outer_vessel_line[2]))])

        if sign(inner_vessel_line[2][1])!=sign(outer_vessel_line[2][1])
            vessel_line=[vcat(inner_vessel_line[1],outer_vessel_line[1]),
                        vcat(inner_vessel_line[2],outer_vessel_line[2])]
        else
            vessel_line=[vcat(inner_vessel_line[1],reverse(outer_vessel_line[1])),
                        vcat(inner_vessel_line[2],reverse(outer_vessel_line[2]))]
        end
        vessel_poly=xy_polygon(vessel_line...)

        # make the divertor domes in the vessel
        cl = Contour.contour(r,z,eqt.profiles_2d[1].psi,eqt.global_quantities.psi_boundary*(1-δψ)+eqt.global_quantities.psi_axis*δψ)
        for line in Contour.lines(cl)
            pr, pz = Contour.coordinates(line)
            if pr[1]!=pr[end]
                vessel_poly=LibGEOS.difference(vessel_poly,xy_polygon(pr,pz))
            end
        end
        
        # cut the top/bottom part of the vessel with the inner_wall_line
        vessel_poly=LibGEOS.intersection(vessel_poly,xy_polygon(inner_wall_line...))
    end

    IMAS.get_radial_build(rb, -1).outline.r = [v[1] for v in LibGEOS.coordinates(vessel_poly)[1]]
    IMAS.get_radial_build(rb, -1).outline.z = [v[2] for v in LibGEOS.coordinates(vessel_poly)[1]]

    IMAS.get_radial_build(rb, -5).outline.r = [v[1] for v in LibGEOS.coordinates(wall_poly)[1]]
    IMAS.get_radial_build(rb, -5).outline.z = [v[2] for v in LibGEOS.coordinates(wall_poly)[1]]

    for (k,layer) in reverse(collect(enumerate(rb.center_stack)))
        if (layer.index<-1) && (layer.index>-5)
            poly=LibGEOS.buffer(xy_polygon(rb.center_stack[k+1].outline.r,rb.center_stack[k+1].outline.z),layer.thickness)
            rb.center_stack[k].outline.r = [v[1] for v in LibGEOS.coordinates(poly)[1]]
            rb.center_stack[k].outline.z = [v[2] for v in LibGEOS.coordinates(poly)[1]]
        end
    end

    IMAS.get_radial_build(rb, 1).outline.r,IMAS.get_radial_build(rb, 1).outline.z = wall_oh(rb)
    IMAS.get_radial_build(rb, 0).outline.r,IMAS.get_radial_build(rb, 0).outline.z = wall_plug(rb)
    rb.center_stack[end].outline.r,rb.center_stack[end].outline.z = wall_cryostat(rb)
    return rb
end

#= == =#
#  OH  #
#= == =#


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
    innerSolenoidRadius, outerSolenoidRadius = (IMAS.get_radial_build(rb, 1).start_radius, IMAS.get_radial_build(rb, 1).end_radius)
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
    R0_TF = sum((IMAS.get_radial_build(rb, -1).start_radius, IMAS.get_radial_build(rb, -1).end_radius)) / 2.0
    Rtf1 = IMAS.get_radial_build(rb, 2).start_radius
    Rtf2 = IMAS.get_radial_build(rb, 2).end_radius
    B0_OH = dd.radial_build.oh_b_field_max
    R_sol1 = IMAS.get_radial_build(rb, 1).start_radius
    R_sol2 = IMAS.get_radial_build(rb, 1).end_radius
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
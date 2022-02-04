import LibGEOS
import Interpolations
import Contour
import LazySets

#= ================== =#
#  init core_profiles  #
#= ================== =#
function init(cp::IMAS.core_profiles; kw...)
    for item in keys(kw)
        IMAS.set_time_array(cp.global_quantities, item, kw[item])
    end
    return cp
end

#= ========== =#
#  init build  #
#= ========== =#
"""
    init(bd::IMAS.build; layers...)

Initialize build IDS based on center stack layers (thicknesses)

NOTE: layer[:].type and layer[:].material follows from naming of layers
*   0 ...gap... : vacuum
*   1 OH: ohmic coil
*   2 TF: toroidal field coil
*   3 shield...: neutron shield
*   4 blanket...: neutron blanket
*   5 wall....: 
*  -1 ...vessel...: 

layer[:].hfs is set depending on if "hfs" or "lfs" appear in the name

layer[:].identifier is created as a hash of then name removing "hfs" or "lfs"
"""
function init(bd::IMAS.build; layers...)
    # empty build IDS
    empty!(bd)
    # assign layers
    resize!(bd.layer, length([layer_name for (layer_name, layer_thickness) in layers if layer_thickness >= 0.0]))
    k = 0
    for (layer_name, layer_thickness) in layers
        if layer_thickness < 0.0
            continue
        end
        k += 1
        bd.layer[k].thickness = layer_thickness
        bd.layer[k].name = replace(String(layer_name), "_" => " ")
        if occursin("gap", lowercase(bd.layer[k].name))
            bd.layer[k].type = 0
            bd.layer[k].material = "vacuum"
        elseif uppercase(bd.layer[k].name) == "OH"
            bd.layer[k].type = 1
        elseif occursin("TF", uppercase(bd.layer[k].name))
            bd.layer[k].type = 2
        elseif occursin("shield", lowercase(bd.layer[k].name))
            bd.layer[k].type = 3
        elseif occursin("blanket", lowercase(bd.layer[k].name))
            bd.layer[k].type = 4
        elseif occursin("wall", lowercase(bd.layer[k].name))
            bd.layer[k].type = 5
        end
        if occursin("hfs", lowercase(bd.layer[k].name))
            bd.layer[k].hfs = 1
        elseif occursin("lfs", lowercase(bd.layer[k].name))
            bd.layer[k].hfs = -1
        else
            bd.layer[k].hfs = 0
        end
        if occursin("vessel", lowercase(bd.layer[k].name))
            bd.layer[k].type = -1
            bd.layer[k].material = "vacuum"
        end
        bd.layer[k].identifier = UInt(hash(replace(replace(lowercase(bd.layer[k].name), "hfs" => ""), "lfs" => "")))
    end
    if ismissing(bd.layer[end],:material) || bd.layer[end].material != "vacuum"
        error("Material of last layer ($(bd.layer[end].name)) must be `vacuum`")
    end

    return bd
end

"""
    init(bd::IMAS.build, eqt::IMAS.equilibrium; is_nuclear_facility=true)

Simple initialization of build IDS based on equilibrium time_slice
"""
function init(bd::IMAS.build,
              eq::IMAS.equilibrium;
              tf_shape_index::Int=3,
              is_nuclear_facility::Bool=true,
              pf_inside_tf::Bool=false,
              pf_outside_tf::Bool=true)
    eqt = eq.time_slice[]
    rmin = eqt.boundary.geometric_axis.r - eqt.boundary.minor_radius
    rmax = eqt.boundary.geometric_axis.r + eqt.boundary.minor_radius

    if is_nuclear_facility
        n_hfs_layers = 6
        gap = (rmax - rmin) / 20.0 # plasma-wall gap
        rmin -= gap
        rmax += gap
        dr = rmin / n_hfs_layers
        init(bd,
            gap_OH=dr * 2.0,
            OH=dr,
            hfs_TF=dr,
            gap_hfs_TF_shield=pf_inside_tf ? 0 : -1,
            hfs_shield=dr / 2.0,
            hfs_blanket=dr,
            hfs_wall=dr / 2.0,
            vacuum_vessel=rmax - rmin,
            lfs_wall=dr / 2.0,
            lfs_blanket=dr * 2,
            lfs_shield=dr / 2.0,
            gap_lfs_TF_shield=dr * (pf_inside_tf ? 4 : -1),
            lfs_TF=dr,
            gap_cryostat=dr * (pf_outside_tf ? 5 : 1))

    else
        n_hfs_layers = 4.5
        gap = (rmax - rmin) / 20.0 # plasma-wall gap
        rmin -= gap
        rmax += gap
        dr = rmin / n_hfs_layers
        init(bd,
            gap_OH=dr * 2.0,
            OH=dr,
            hfs_TF=dr,
            gap_hfs_TF_wall=0.0,
            hfs_wall=dr / 2.0,
            vacuum_vessel=rmax - rmin,
            lfs_wall=dr / 2.0,
            gap_lfs_TF_wall=dr * 3,
            lfs_TF=dr,
            gap_cryostat=2 * dr)
    end

    # TF coils
    bd.tf.coils_n = 16

    # cross-section outlines
    build_cx(bd, eqt, tf_shape_index)

    return bd
end

function wall_miller_conformal(bd, layer_type, elongation, triangularity; n_points=101)
    if layer_type == -1
        Rstart = IMAS.get_build(bd, type=layer_type).start_radius
        Rend = IMAS.get_build(bd, type=layer_type).end_radius
        line = miller_Rstart_Rend(Rstart, Rend, elongation, triangularity; n_points)        
        return line, line
    else
        Rstart_lfs = IMAS.get_build(bd, type=layer_type, hfs=-1).start_radius
        Rend_lfs = IMAS.get_build(bd, type=layer_type, hfs=-1).end_radius
        Rstart_hfs = IMAS.get_build(bd, type=layer_type, hfs=1).start_radius
        Rend_hfs = IMAS.get_build(bd, type=layer_type, hfs=1).end_radius
        inner_line = miller_Rstart_Rend(Rend_hfs, Rstart_lfs, elongation, triangularity; n_points)
        outer_line = miller_Rstart_Rend(Rstart_hfs, Rend_lfs, elongation, triangularity; n_points)
        return inner_line, outer_line
    end
end

"""
    build_cx(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice)
Translates 1D build to 2D cross-sections
"""

function build_cx(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice, tf_shape_index::Int)
    # Inner radii of the vacuum vessel
    R_hfs_vessel = IMAS.get_build(bd, type=-1).start_radius
    R_lfs_vessel = IMAS.get_build(bd, type=-1).end_radius
    
    ψb = IMAS.find_psi_boundary(eqt)
    ψa = eqt.profiles_1d.psi[1]
    # Vessel as buffered convex-hull polygon of LCFS and strike points
    r95, z95, _ = IMAS.flux_surface(eqt, (ψb - ψa) * 0.90 + ψa, true)
    Z0 = eqt.global_quantities.magnetic_axis.z
    rlcfs, zlcfs, _ = IMAS.flux_surface(eqt, ψb, true)
    theta = range(0.0, 2 * pi, length=101)
    private_extrema = []
    private = IMAS.flux_surface(eqt, ψb, false)
    a = 0
    for (pr, pz) in private
        if sign(pz[1] - Z0) != sign(pz[end] - Z0)
            # open flux surface does not encicle the plasma
            continue
        elseif minimum_distance_two_shapes(pr, pz, rlcfs, zlcfs) > (maximum(zlcfs) - minimum(zlcfs)) / 20
            # secondary Xpoint far away
            continue
        elseif (sum(pz) - Z0) < 0
            # lower private region
            index = argmax(pz)
            a = minimum(z95) - minimum(zlcfs)
            a = min(a, pz[index] - minimum(pz))
        else
            # upper private region
            index = argmin(pz)
            a = maximum(zlcfs)-maximum(z95)
            a = min(a, maximum(pz) - pz[index])
        end
        Rx = pr[index]
        Zx = pz[index]
        append!(private_extrema, IMAS.intersection(a .* cos.(theta) .+ Rx, a .* sin.(theta) .+ Zx, pr, pz))
    end
    h = [[r,z] for (r,z) in vcat(collect(zip(rlcfs,zlcfs)),private_extrema)]
    hull = LazySets.convex_hull(h)
    R = [r for (r,z) in hull]
    R[R .< R_hfs_vessel] .= R_hfs_vessel
    R[R .> R_lfs_vessel] .= R_lfs_vessel
    Z = [z for (r,z) in hull]
    hull_poly = xy_polygon(R, Z)
    vessel_poly = LibGEOS.buffer(hull_poly, ((R_lfs_vessel - R_hfs_vessel) - (maximum(rlcfs) - minimum(rlcfs))) / 2.0 )

    # make the divertor domes in the vessel
    δψ = 0.05
    r = range(eqt.profiles_2d[1].grid.dim1[1], eqt.profiles_2d[1].grid.dim1[end], length=length(eqt.profiles_2d[1].grid.dim1))
    z = range(eqt.profiles_2d[1].grid.dim2[1], eqt.profiles_2d[1].grid.dim2[end], length=length(eqt.profiles_2d[1].grid.dim2))
    cl = Contour.contour(r, z, eqt.profiles_2d[1].psi, eqt.profiles_1d.psi[end] * (1 - δψ) + eqt.profiles_1d.psi[1] * δψ)
    for line in Contour.lines(cl)
        pr, pz = Contour.coordinates(line)
        if pr[1] != pr[end]
            pz[1] = pz[1] * 2
            pz[end] = pz[end] * 2
            vessel_poly = LibGEOS.difference(vessel_poly, xy_polygon(pr, pz))
        end
    end

    # vacuum vessel
    IMAS.get_build(bd, type=-1).outline.r = [v[1] for v in LibGEOS.coordinates(vessel_poly)[1]]
    IMAS.get_build(bd, type=-1).outline.z = [v[2] for v in LibGEOS.coordinates(vessel_poly)[1]]

    # all layers between vessel and OH
    vessel_to_oh = []
    valid = false
    for (k, layer) in reverse(collect(enumerate(bd.layer)))
        # stop once you see the OH
        if layer.type == 1
            valid = false
            break
        end
        if valid
            push!(vessel_to_oh, k)
        end
        # valid starting from the vessel
        if layer.type == -1
            valid = true
        end
    end
    shape_set = false
    for (n, k) in enumerate(vessel_to_oh)
        # layer that preceeds the TF (or shield) sets the TF (and shield) shape
        if (!shape_set) && (n<length(vessel_to_oh)) && (bd.layer[vessel_to_oh[n+1]].type in [2, 3])
            FUSE.optimize_shape(bd, k, tf_shape_index)
            shape_set = true
        # everything else is conformal convex hull
        else
            FUSE.optimize_shape(bd, k, -2)
        end
    end

    # plug
    L = 0
    R = IMAS.get_build(bd, type=1).start_radius
    D = minimum(IMAS.get_build(bd, type=2, hfs=1).outline.z)
    U = maximum(IMAS.get_build(bd, type=2, hfs=1).outline.z)
    bd.layer[1].outline.r, bd.layer[1].outline.z = rectangle_shape(L, R, D, U)

    # oh
    L = IMAS.get_build(bd, type=1).start_radius
    R = IMAS.get_build(bd, type=1).end_radius
    D = minimum(IMAS.get_build(bd, type=2, hfs=1).outline.z)
    U = maximum(IMAS.get_build(bd, type=2, hfs=1).outline.z)
    bd.layer[2].outline.r, bd.layer[2].outline.z = rectangle_shape(L, R, D, U)

    # cryostat
    L = 0
    R = bd.layer[end].end_radius
    D = minimum(IMAS.get_build(bd, type=2, hfs=1).outline.z) - bd.layer[end].thickness
    U = maximum(IMAS.get_build(bd, type=2, hfs=1).outline.z) + bd.layer[end].thickness
    bd.layer[end].outline.r, bd.layer[end].outline.z = rectangle_shape(L, R, D, U)

    # set the toroidal thickness of the TF coils based on the innermost radius and the number of coils
    bd.tf.thickness = 2 * π * IMAS.get_build(bd, type=2, hfs=1).start_radius / bd.tf.coils_n
    return bd
end


function optimize_shape(bd::IMAS.build, layer_index::Int, tf_shape_index::Int)
    # properties of current layer
    layer = bd.layer[layer_index]
    id = bd.layer[layer_index].identifier
    r_start = layer.start_radius
    r_end = IMAS.get_build(bd, identifier=layer.identifier, hfs=-1).end_radius
    hfs_thickness = layer.thickness
    lfs_thickness = IMAS.get_build(bd, identifier=id, hfs=-1).thickness
    target_minimum_distance = (hfs_thickness + lfs_thickness) / 2.0

    # obstruction
    oR = bd.layer[layer_index+1].outline.r
    oZ = bd.layer[layer_index+1].outline.z

    if ismissing(layer, :shape)
        layer.shape = tf_shape_index
    end

    # handle offset and offset & convex-hull
    if layer.shape in [-1, -2]
        poly = LibGEOS.buffer(xy_polygon(oR, oZ), (hfs_thickness + lfs_thickness) / 2.0)
        layer.outline.r = [v[1] .+ (lfs_thickness .- hfs_thickness) / 2.0 for v in LibGEOS.coordinates(poly)[1]]
        layer.outline.z = [v[2] for v in LibGEOS.coordinates(poly)[1]]
        if layer.shape == -2
            h = [[r,z] for (r,z) in collect(zip(layer.outline.r,layer.outline.z))]
            hull =LazySets.convex_hull(h)
            layer.outline.r = vcat([r for (r,z) in hull],hull[1][1])
            layer.outline.z = vcat([z for (r,z) in hull],hull[1][2])
        end
    # handle shapes
    else
        up_down_symmetric = false
        if abs(sum(oZ)/sum(abs.(oZ))) < 1E-2
            up_down_symmetric = true
        end

        if up_down_symmetric
            layer.shape = mod(layer.shape, 100)
        else
            layer.shape = mod(layer.shape, 100) + 100
        end
        
        func = shape_function(layer.shape)
        if ismissing(layer, :shape_parameters)
            layer.shape_parameters = init_shape_parameters(layer.shape, oR, oZ, r_start, r_end, target_minimum_distance)
        end
        layer.outline.r, layer.outline.z = func(r_start, r_end, layer.shape_parameters...)
        layer.shape_parameters = optimize_shape(oR, oZ, target_minimum_distance, func, r_start, r_end, layer.shape_parameters)
        layer.outline.r, layer.outline.z = func(r_start, r_end, layer.shape_parameters...)
    end
end

#= ================ =#
#  flux-swing actor #
#= ================ =#

mutable struct FluxSwingActor <: AbstractActor
    bd::IMAS.build
    eqt::IMAS.equilibrium__time_slice
    cp::IMAS.core_profiles
end

function FluxSwingActor(dd::IMAS.dd)
    return FluxSwingActor(dd.build, dd.equilibrium, dd.core_profiles)
end

function FluxSwingActor(bd::IMAS.build, eq::IMAS.equilibrium, cp::IMAS.core_profiles)
    return FluxSwingActor(bd, eq.time_slice[], cp)
end

# step
function step(flxactor::FluxSwingActor)
    rampup_flux_requirements(flxactor.bd, flxactor.eqt, flxactor.cp)
    flattop_flux_requirements(flxactor.bd, flxactor.eqt, flxactor.cp)
    pf_flux_requirements(flxactor.bd, flxactor.eqt)
end

"""
    rampup_flux_requirements(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice, cp::IMAS.core_profiles)

Estimate OH flux requirement during rampup

NOTES:
* Equations from GASC (Stambaugh FST 2011)
* eqt is supposed to be the equilibrium right at the end of the rampup phase, beginning of flattop
* core_profiles is only used to get core_profiles.global_quantities.ejima
"""
function rampup_flux_requirements(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice, cp::IMAS.core_profiles)
    # from IMAS dd to local variables
    majorRadius = eqt.boundary.geometric_axis.r
    minorRadius = eqt.boundary.minor_radius
    elongation = eqt.boundary.elongation
    plasmaCurrent = eqt.global_quantities.ip / 1E6 # in [MA]
    li = eqt.global_quantities.li_3 # what li ?
    ejima = @ddtime cp.global_quantities.ejima

    # ============================= #
    # evaluate plasma inductance
    plasmaInductanceInternal = 0.4 * 0.5 * pi * majorRadius * li
    plasmaInductanceExternal = 0.4 * pi * majorRadius * (log(8.0 * majorRadius / minorRadius / sqrt(elongation)) - 2.0)
    plasmaInductanceTotal = plasmaInductanceInternal + plasmaInductanceExternal

    # estimate rampup flux requirement
    rampUpFlux = (ejima * 0.4 * pi * majorRadius + plasmaInductanceTotal) * plasmaCurrent

    # ============================= #
    bd.flux_swing_requirements.rampup = abs(rampUpFlux)
end

"""
    flattop_flux_requirements(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice, cp::IMAS.core_profiles)

Estimate OH flux requirement during flattop 

NOTES:
* this is a dummy function right now!, we simply take 1/2 of the rampup
"""
function flattop_flux_requirements(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice, cp::IMAS.core_profiles)
    # from IMAS dd to local variables

    # ============================= #
    #plasmaResistivity = calc_plasmaResitivity(IN['plasma parameters']['Ti0'],neLinAvg,effectiveZ,Tratio,St,Sn,1.0,1./aspectRatio)
    #flattopFluxConsumption = 1.e6*plasmaResistivity * inductiveFraction * plasmaCurrent * flattopDuration

    flattopFluxConsumption = 0.5 * bd.flux_swing_requirements.rampup

    # ============================= #
    bd.flux_swing_requirements.flattop = flattopFluxConsumption
end

"""
    pf_flux_requirements(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice)

Estimate vertical field from PF coils and its contribution to flux swing

NOTES:
* Equations from GASC (Stambaugh FST 2011)
* eqt is supposed to be the equilibrium right at the end of the rampup phase, beginning of flattop
"""
function pf_flux_requirements(bd::IMAS.build, eqt::IMAS.equilibrium__time_slice)
    # from IMAS dd to local variables
    majorRadius = eqt.boundary.geometric_axis.r
    minorRadius = eqt.boundary.minor_radius
    elongation = eqt.boundary.elongation
    plasmaCurrent = eqt.global_quantities.ip / 1E6 # in [MA]
    betaP = eqt.global_quantities.beta_pol
    li = eqt.global_quantities.li_3 # what li ?

    # ============================= #
    # estimate vertical field and its contribution to flux swing
    verticalFieldAtCenter = 0.1 * plasmaCurrent / majorRadius * (log(8.0 * majorRadius / (minorRadius * sqrt(elongation))) - 1.5 + betaP + 0.5 * li)
    fluxFromVerticalField = 0.8 * verticalFieldAtCenter * pi * (majorRadius^2 - (majorRadius - minorRadius)^2)

    # ============================= #
    bd.flux_swing_requirements.pf = - abs(fluxFromVerticalField)
end

"""
    oh_requirements(bd::IMAS.build, double_swing::Bool=true)

Evaluate OH current density and B_field required for rampup and flattop

NOTES:
* Equations from GASC (Stambaugh FST 2011)
* Also relevant: `Engineering design solutions of flux swing with structural requirements for ohmic heating solenoids` Smith, R. A. September 30, 1977
"""
function oh_requirements(bd::IMAS.build, double_swing::Bool=true)
    innerSolenoidRadius, outerSolenoidRadius = (IMAS.get_build(bd, type=1).start_radius, IMAS.get_build(bd, type=1).end_radius)
    totalOhFluxReq = bd.flux_swing_requirements.rampup.total + bd.flux_swing_requirements.flattop + bd.flux_swing_requirements.pf

    # ============================= #

    # Calculate magnetic field at solenoid bore required to match flux swing request
    RiRoFactor = innerSolenoidRadius / outerSolenoidRadius
    magneticFieldSolenoidBore = 3.0 * totalOhFluxReq / pi / outerSolenoidRadius^2 / (RiRoFactor^2 + RiRoFactor + 1.0) / (double_swing ? 2 : 1)
    currentDensityOH = magneticFieldSolenoidBore / (0.4 * pi * outerSolenoidRadius*(1-innerSolenoidRadius/outerSolenoidRadius))

    # ============================= #

    # minimum requirements for OH
    bd.oh.required.b_field = magneticFieldSolenoidBore
    bd.oh.required.j = currentDensityOH
end

#= ======== =#
#  Stresses  #
#= ======== =#

function stress_calculations(dd::IMAS.dd)
    error("not completed yet")
    bd = dd.build
    B0_TF = dd.build.tf_b_field_max
    R0_TF = sum((IMAS.get_build(bd, type=-1).start_radius, IMAS.get_build(bd, type=-1).end_radius)) / 2.0
    Rtf1 = IMAS.get_build(bd, type=2).start_radius
    Rtf2 = IMAS.get_build(bd, type=2).end_radius
    B0_OH = dd.build.oh_b_field_max
    R_sol1 = IMAS.get_build(bd, type=1).start_radius
    R_sol2 = IMAS.get_build(bd, type=1).end_radius
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
#= ============== =#
#  OH TF stresses  #
#= ============== =#
Base.@kwdef mutable struct FUSEparameters__ActorStresses{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    do_plot::Entry{Bool} = Entry{Bool}("-", "Plot"; default=false)
    n_points::Entry{Int} = Entry{Int}("-", "Number of grid points"; default=5)
end

mutable struct ActorStresses{D,P} <: ReactorAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorStresses{P}
    function ActorStresses(dd::IMAS.dd{D}, par::FUSEparameters__ActorStresses{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorStresses)
        par = par(kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorStresses(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates mechanical stresses on the center stack

!!! note 
    Stores data in `dd.solid_mechanics`
"""
function ActorStresses(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorStresses(dd, act.ActorStresses; kw...)
    step(actor)
    finalize(actor)
    if actor.par.do_plot
        display(plot(actor.dd.solid_mechanics.center_stack.stress))
    end
    return actor
end

function _step(actor::ActorStresses)
    dd = actor.dd
    par = actor.par

    eq = dd.equilibrium
    bd = dd.build
    sm = dd.solid_mechanics

    plasma = IMAS.get_build_layer(bd.layer, type=_plasma_)
    R0 = (plasma.end_radius + plasma.start_radius) / 2.0
    B0 = maximum(abs, eq.vacuum_toroidal_field.b0)

    R_tf_in = IMAS.get_build_layer(bd.layer, type=_tf_, fs=_hfs_).start_radius
    R_tf_out = IMAS.get_build_layer(bd.layer, type=_tf_, fs=_hfs_).end_radius
    
    Bz_oh = bd.oh.max_b_field
    
    R_oh_in = IMAS.get_build_layer(bd.layer, type=_oh_).start_radius
    R_oh_out = IMAS.get_build_layer(bd.layer, type=_oh_).end_radius
    
    f_struct_tf = bd.tf.technology.fraction_steel
    f_struct_oh = bd.oh.technology.fraction_steel

    bucked = sm.center_stack.bucked == 1
    noslip = sm.center_stack.noslip == 1
    plug = sm.center_stack.plug == 1
    
    for oh_on in (true, false)
        solve_1D_solid_mechanics!(
            sm.center_stack,
            R0,
            B0,
            R_tf_in,
            R_tf_out,
            oh_on ? Bz_oh : 0.0,
            R_oh_in,
            R_oh_out;
            bucked=bucked,
            noslip=noslip,
            plug=plug,
            f_struct_tf=f_struct_tf,
            f_struct_oh=f_struct_oh,
            f_struct_pl=1.0,
            em_tf=sm.center_stack.properties.young_modulus.tf,
            gam_tf=sm.center_stack.properties.poisson_ratio.tf,
            em_oh=sm.center_stack.properties.young_modulus.oh,
            gam_oh=sm.center_stack.properties.poisson_ratio.oh,
            em_pl=getproperty(sm.center_stack.properties.young_modulus, :pl, NaN),
            gam_pl=getproperty(sm.center_stack.properties.poisson_ratio, :pl, NaN),
            n_points=par.n_points,
            empty_smcs=oh_on,
            verbose=false
        )
    end

    return actor
end

@recipe function plot_ActorStresses(actor::ActorStresses)
    @series begin
        actor.dd.solid_mechanics.center_stack.stress
    end
end

"""
    solve_1D_solid_mechanics!(
        smcs::IMAS.solid_mechanics__center_stack,
        R0::T,                                 # : (float) major radius at center of TF bore, meters
        B0::T,                                 # : (float) toroidal field at R0, Tesla
        R_tf_in::T,                            # : (float) major radius of inboard edge of TF coil core legs, meters
        R_tf_out::T,                           # : (float) major radius of outboard edge of TF coil core legs, meters
        Bz_oh::T,                              # : (float) axial field in solenoid bore, Tesla
        R_oh_in::T,                            # : (float) major radius of inboard edge of OH coil, meters
        R_oh_out::T;                           # : (float) major radius of outboard edge of OH coil, meters
        axial_stress_tf_avg::T=NaN,            # : (float) average axial stress in TF coil core legs, Pa (if nothing, use constant fraction of hoop stress)
        axial_stress_oh_avg::T=NaN,            # : (float) average axial stress in OH coil, Pa (if nothing, use constant fraction of hoop stress)
        bucked::Bool,                          # : (bool), flag for bucked boundary conditions between TF and OH (and center plug, if present)
        noslip::Bool,                          # : (bool), flag for no slip conditions between TF and OH (and center plug, if present)
        plug::Bool,                            # : (bool), flag for center plug
        f_struct_tf::T,                        # : (float), fraction of TF coil that is structural material
        f_struct_oh::T,                        # : (float), fraction of OH coil that is structural material
        f_struct_pl::T,                        # : (float), fraction of plug that is structural material
        em_tf::T,                              # : (float), modulus of elasticity for TF coil, Pa (default is stainless steel)
        gam_tf::T,                             # : (float), Poisson`s ratio for TF coil, (default is stainless steel)
        em_oh::T,                              # : (float), modulus of elasticity for OH coil, Pa (default is stainless steel)
        gam_oh::T,                             # : (float), Poisson`s ratio for OH coil, (default is stainless steel)
        em_pl::T,                              # : (float), modulus of elasticity for center plug, Pa (default is stainless steel)
        gam_pl::T,                             # : (float), Poisson`s ratio for center plug, (default is stainless steel)
        f_tf_sash::T=0.873,                    # : (float), conversion factor from hoop stress to axial stress for TF coil
        f_oh_sash::T=0.37337,                  # : (float), conversion factor from hoop stress to axial stress for OH coil
        n_points::Int=21,                      # : (int), number of radial points
        empty_smcs::Bool=true,                 # : (bool), flag to empty the smcs structure (useful to identify worst case in a series of scenarios)
        verbose::Bool=false) where (T<:Real)   # : (bool), flag for verbose output to terminal

Uses Leuer 1D solid mechanics equations to solve for radial and hoop stresses in TF coil, OH coil, and center plug.
Based on derivations in Engineering Physics Note "EPNjal17dec17_gasc_pt5_tf_oh_plug_buck" by Jim Leuer (Dec. 17, 2017)

Returns radial, hoop, axial, and Von Mises stresses for TF, OH, and plug (Pascals)
(optional) radial profiles of radial, hoop, axial, and Von Mises stresses for TF, OH, and plug (Pascals)

The tokamak radial buid is :

|| plug or void (0 < r < R1) || coil 1 (R1 < r < R2) || coil 2 (R3 < r < R4) || ----> plasma center (r = R0)
"""
function solve_1D_solid_mechanics!(
    smcs::IMAS.solid_mechanics__center_stack,
    R0::T,                                 # : (float) major radius at center of TF bore, meters
    B0::T,                                 # : (float) toroidal field at R0, Tesla
    R_tf_in::T,                            # : (float) major radius of inboard edge of TF coil core legs, meters
    R_tf_out::T,                           # : (float) major radius of outboard edge of TF coil core legs, meters
    Bz_oh::T,                              # : (float) axial field in solenoid bore, Tesla
    R_oh_in::T,                            # : (float) major radius of inboard edge of OH coil, meters
    R_oh_out::T;                           # : (float) major radius of outboard edge of OH coil, meters
    axial_stress_tf_avg::T=NaN,            # : (float) average axial stress in TF coil core legs, Pa (if nothing, use constant fraction of hoop stress)
    axial_stress_oh_avg::T=NaN,            # : (float) average axial stress in OH coil, Pa (if nothing, use constant fraction of hoop stress)
    bucked::Bool,                          # : (bool), flag for bucked boundary conditions between TF and OH (and center plug, if present)
    noslip::Bool,                          # : (bool), flag for no slip conditions between TF and OH (and center plug, if present)
    plug::Bool,                            # : (bool), flag for center plug
    f_struct_tf::T,                        # : (float), fraction of TF coil that is structural material
    f_struct_oh::T,                        # : (float), fraction of OH coil that is structural material
    f_struct_pl::T,                        # : (float), fraction of plug that is structural material
    em_tf::T,                              # : (float), modulus of elasticity for TF coil, Pa (default is stainless steel)
    gam_tf::T,                             # : (float), Poisson`s ratio for TF coil, (default is stainless steel)
    em_oh::T,                              # : (float), modulus of elasticity for OH coil, Pa (default is stainless steel)
    gam_oh::T,                             # : (float), Poisson`s ratio for OH coil, (default is stainless steel)
    em_pl::T,                              # : (float), modulus of elasticity for center plug, Pa (default is stainless steel)
    gam_pl::T,                             # : (float), Poisson`s ratio for center plug, (default is stainless steel)
    f_tf_sash::T=0.873,                    # : (float), conversion factor from hoop stress to axial stress for TF coil
    f_oh_sash::T=0.37337,                  # : (float), conversion factor from hoop stress to axial stress for OH coil
    n_points::Int=21,                      # : (int), number of radial points
    empty_smcs::Bool=true,                 # : (bool), flag to empty the smcs structure (useful to identify worst case in a series of scenarios)
    verbose::Bool=false) where (T<:Real)   # : (bool), flag for verbose output to terminal

    if empty_smcs
        # empty smcs but always retain the materials' properties
        for field in keys(smcs)
            if field != :properties
                empty!(smcs, field)
            end
        end
    end

    tp = typeof(promote(R0, B0, R_tf_in, R_tf_out, Bz_oh, R_oh_in, R_oh_out)[1])

    if verbose
        println("solve_1D_solid_mechanics:")
        println("- R0 = ", R0)
        println("- B0 = ", B0)
        println("- R_tf_in = ", R_tf_in)
        println("- R_tf_out = ", R_tf_out)
        println("- Bz_oh = ", Bz_oh)
        println("- R_oh_in = ", R_oh_in)
        println("- R_oh_out = ", R_oh_out)
        println("- axial_stress_tf_avg = ", axial_stress_tf_avg)
        println("- axial_stress_oh_avg = ", axial_stress_oh_avg)
        println("- bucked = ", bucked)
        println("- noslip = ", noslip)
        println("- plug = ", plug)
        println("- f_struct_tf = ", f_struct_tf)
        println("- f_struct_oh = ", f_struct_oh)
        println("- f_struct_pl = ", f_struct_pl)
    end

    # define structural constants
    embar_tf = em_tf / (1 - gam_tf^2)
    embar_oh = em_oh / (1 - gam_oh^2)
    embar_pl = em_pl / (1 - gam_pl^2)

    # define forcing constraints on TF and OH coils
    C_tf = 1.0 / embar_tf * 2 * (B0 * R0)^2 / (constants.μ_0 * (R_tf_out^2 - R_tf_in^2)^2)
    C_oh = -1.0 / embar_oh * Bz_oh^2 / (constants.μ_0 * (R_oh_out - R_oh_in)^2)

    # calculate centerlines, check radial build inputs for consistency 
    cl_tf = 0.5 * (R_tf_out + R_tf_in)
    cl_oh = 0.5 * (R_oh_out + R_oh_in)

    # determine ordering of radial build
    if verbose
        if cl_oh < cl_tf
            println("* order = OH-TF")
        else
            println("* order = TF-OH")
        end
    end

    # check radial build inputs for consistency 
    if R_tf_out < R_tf_in
        error("R_tf_out ($R_tf_out) < R_tf_in ($R_tf_in)")
    end
    if R_oh_out < R_oh_in
        error("R_oh_out ($R_oh_out) < R_oh_in ($R_oh_in)")
    end
    if (cl_oh > cl_tf) && ((R_tf_out - R_oh_in * (1 + 1e-3)) > 0)
        error("* TF and OH are overlapping: R_oh_in ($R_oh_in) < R_tf_out ($R_tf_out)")
    end
    if (cl_oh < cl_tf) && ((R_oh_out - R_tf_in * (1 + 1e-3)) > 0)
        error("* TF and OH are overlapping: R_tf_in ($R_tf_in) < R_oh_out ($R_oh_out)")
    end
    if bucked && (cl_oh < cl_tf) && (abs.(R_oh_out - R_tf_in) > R_tf_in * 1e-3)
        error("* OH is bucked against TF, but R_oh_out ($R_oh_out) != R_tf_in ($R_tf_in)")
    end

    # define radial functions in u/r and du/dr for TF and OH
    # u_tf/r   = C_tf * f_tf(r) + A_tf + B_tf/r^2
    # du_tf/dr = C_tf * g_tf(r) + A_tf - B_tf/r^2
    # u_oh/r   = C_oh * f_oh(r) + A_oh + B_oh/r^2
    # du_oh/dr = C_oh * g_oh(r) + A_oh - B_oh/r^2
    function f_tf(r)
        logr = log(r)
        if !isfinite(logr)
            logr = -100.0
        end
        return 1.0 / 8.0 * r^2 - R_tf_in^2 / 2.0 * (logr - 0.5)
    end

    function g_tf(r)
        logr = log(r)
        if !isfinite(logr)
            logr = -100.0
        end
        return 3.0 / 8.0 * r^2 - R_tf_in^2 / 2.0 * (logr + 0.5)
    end

    function f_oh(r)
        return 1.0 / 3.0 * R_oh_out * r - 1.0 / 8.0 * r^2
    end
    function g_oh(r)
        return 2.0 / 3.0 * R_oh_out * r - 3.0 / 8.0 * r^2
    end

    # define linear system of equations based on boundary conditions
    # M * X = Y

    if !bucked
        ## free standing OH and TF coils
        # radial stress = 0 at ALL coil edges
        if plug
            error("TF and OH must be bucked to have a central plug")
        end
        if verbose
            println("* Free standing coils")
        end

        M = zeros(tp, 4, 4)
        M[1, :] = [1 + gam_tf, (gam_tf - 1) / R_tf_in^2, 0.0, 0.0]
        M[2, :] = [1 + gam_tf, (gam_tf - 1) / R_tf_out^2, 0.0, 0.0]
        M[3, :] = [0.0, 0.0, 1 + gam_oh, (gam_oh - 1) / R_oh_in^2]
        M[4, :] = [0.0, 0.0, 1 + gam_oh, (gam_oh - 1) / R_oh_out^2]
        Y = [
            -C_tf * (g_tf(R_tf_in) + gam_tf * f_tf(R_tf_in)),
            -C_tf * (g_tf(R_tf_out) + gam_tf * f_tf(R_tf_out)),
            -C_oh * (g_oh(R_oh_in) + gam_oh * f_oh(R_oh_in)),
            -C_oh * (g_oh(R_oh_out) + gam_oh * f_oh(R_oh_out)),
        ]
        A_tf, B_tf, A_oh, B_oh = M \ Y
        A_pl = 0.0

    elseif !plug && (cl_oh < cl_tf)
        ## bucked OH and TF only (no plug)
        # order must be "OH-TF"
        # radial stress = 0 at R_oh_in and R_tf_out
        # radial stresses and displacementes are equal at interface R_int ( == R_oh_out = R_tf_in )
        if verbose
            println("* OH bucked against TF only")
        end
        R_int = R_oh_out
        M = zeros(tp, 4, 4)
        M[1, :] = [1 + gam_tf, (gam_tf - 1) / R_tf_out^2, 0.0, 0.0]
        M[2, :] = [0.0, 0.0, 1 + gam_oh, (gam_oh - 1) / R_oh_in^2]
        M[3, :] = [1.0, 1.0 / R_int^2, -1.0, -1.0 / R_int^2]
        M[4, :] = [embar_tf * (1 + gam_tf), embar_tf * (gam_tf - 1) / R_int^2, -embar_oh * (1 + gam_oh), -embar_oh * (gam_oh - 1) / R_int^2]

        Y = [
            -C_tf * (g_tf(R_tf_out) + gam_tf * f_tf(R_tf_out)),
            -C_oh * (g_oh(R_oh_in) + gam_oh * f_oh(R_oh_in)),
            C_oh * f_oh(R_int) - C_tf * f_tf(R_int),
            embar_oh * C_oh * (g_oh(R_int) + gam_oh * f_oh(R_int)) - embar_tf * C_tf * (g_tf(R_int) + gam_tf * f_tf(R_int)),
        ]
        A_tf, B_tf, A_oh, B_oh = M \ Y
        A_pl = 0.0

    elseif plug && (cl_oh < cl_tf)
        ## bucked plug, OH, and TF
        # order is "OH-TF"
        # radial stress = 0 at R_tf_out
        # radial stresses and displacements are equal at plug-OH interface R_pl ( == R_oh_in) 
        #  and OH-TF interface R_int ( == R_oh_out = R_tf_in )
        if verbose
            println("* OH bucked against TF and plug")
        end
        R_int = R_oh_out
        R_pl = R_oh_in
        M = zeros(tp, 5, 5)
        M[1, :] = [1 + gam_tf, (gam_tf - 1) / R_tf_out^2, 0.0, 0.0, 0.0]
        M[2, :] = [1.0, 1.0 / R_int^2, -1.0, -1.0 / R_int^2, 0.0]
        M[3, :] = [embar_tf * (1 + gam_tf), embar_tf * (gam_tf - 1) / R_int^2, -embar_oh * (1 + gam_oh), -embar_oh * (gam_oh - 1) / R_int^2, 0.0]
        M[4, :] = [0.0, 0.0, 1.0, 1.0 / R_pl^2, -1.0]
        M[5, :] = [0.0, 0.0, embar_oh * (1 + gam_oh), embar_oh * (gam_oh - 1) / R_pl^2, -embar_pl * (1 + gam_pl)]
        Y = [
            -C_tf * (g_tf(R_tf_out) + gam_tf * f_tf(R_tf_out)),
            C_oh * f_oh(R_int) - C_tf * f_tf(R_int),
            embar_oh * C_oh * (g_oh(R_int) + gam_oh * f_oh(R_int)) - embar_tf * C_tf * (g_tf(R_int) + gam_tf * f_tf(R_int)),
            -C_oh * f_oh(R_pl),
            -embar_oh * C_oh * (g_oh(R_pl) + gam_oh * f_oh(R_pl)),
        ]
        A_tf, B_tf, A_oh, B_oh, A_pl = M \ Y

    elseif plug && (cl_oh > cl_tf)
        ## bucked TF aginst plug, OH is free standing
        # order is "TF-OH"
        # radial stress = 0 at R_tf_out
        # radial stress = 0 at R_oh_in and R_oh_out
        # radial stresses and displacements are equal at plug-TF interface R_pl ( == R_tf_in)
        if verbose
            println("* TF bucked against plug")
        end
        R_pl = R_tf_in
        M = zeros(tp, 5, 5)
        M[1, :] = [1 + gam_tf, (gam_tf - 1) / R_tf_out^2, 0.0, 0.0, 0.0]
        M[2, :] = [0.0, 0.0, 1 + gam_oh, (gam_oh - 1) / R_oh_in^2, 0.0]
        M[3, :] = [0.0, 0.0, 1 + gam_oh, (gam_oh - 1) / R_oh_out^2, 0.0]
        M[4, :] = [1.0, 1.0 / R_pl^2, 0.0, 0.0, -1.0]
        M[5, :] = [embar_tf * (1 + gam_tf), embar_tf * (gam_tf - 1) / R_pl^2, 0.0, 0.0, -embar_pl * (1 + gam_pl)]
        Y = [
            -C_tf * (g_tf(R_tf_out) + gam_tf * f_tf(R_tf_out)),
            -C_oh * (g_oh(R_oh_in) + gam_oh * f_oh(R_oh_in)),
            -C_oh * (g_oh(R_oh_out) + gam_oh * f_oh(R_oh_out)),
            -C_tf * f_tf(R_pl),
            -embar_tf * C_tf * (g_tf(R_pl) + gam_tf * f_tf(R_pl)),
        ]
        A_tf, B_tf, A_oh, B_oh, A_pl = M \ Y

    else
        error("radial build inputs do not match known boundary conditions!!")
    end

    ## use integration constants to define functions for displacement and radial, hoop, and Vom Mises stresses

    function u_tf(r)
        return @. r * (C_tf * f_tf(r) + A_tf + B_tf / r^2)
    end

    function dudr_tf(r)
        return C_tf * g_tf(r) + A_tf - B_tf / r^2
    end

    function u_oh(r)
        return r * (C_oh * f_oh(r) + A_oh + B_oh / r^2)
    end

    function dudr_oh(r)
        return C_oh * g_oh(r) + A_oh - B_oh / r^2
    end

    function u_pl(r)
        return r * A_pl
    end

    function dudr_pl(r)
        return A_pl
    end

    function sr(r, em, gam, u, dudr)
        return em / (1 - gam^2) * (dudr + u / r * gam)
    end

    function sh(r, em, gam, u, dudr)
        return em / (1 - gam^2) * (u / r + dudr * gam)
    end

    function svm(sr, sh, sa)
        return sqrt(((sh - sa)^2 + (sa - sr)^2 + (sr - sh)^2) / 2.0)
    end

    ## calculate radial and hoop stresses
    # also estimate axial stresses if not given & modify due to noslip condition
    # default axial scalings taken from Bending free formula & GASC
    r_oh = LinRange(R_oh_in, R_oh_out, n_points)
    r_tf = LinRange(R_tf_in, R_tf_out, n_points)

    displacement_oh = u_oh.(r_oh)
    displacement_tf = u_tf.(r_tf)
    ddiplacementdr_oh = dudr_oh.(r_oh)
    ddiplacementdr_tf = dudr_tf.(r_tf)

    radial_stress_oh = sr.(r_oh, em_oh, gam_oh, displacement_oh, ddiplacementdr_oh)
    radial_stress_tf = sr.(r_tf, em_tf, gam_tf, displacement_tf, ddiplacementdr_tf)

    hoop_stress_oh = sh.(r_oh, em_oh, gam_oh, displacement_oh, ddiplacementdr_oh)
    hoop_stress_tf = sh.(r_tf, em_tf, gam_tf, displacement_tf, ddiplacementdr_tf)

    if isnan(axial_stress_tf_avg)
        hoop_stress_tf_avg = sum(hoop_stress_tf) / length(hoop_stress_tf)
        axial_stress_tf_avg = -f_tf_sash * hoop_stress_tf_avg
    end
    if isnan(axial_stress_oh_avg)
        hoop_stress_oh_avg = sum(hoop_stress_oh) / length(hoop_stress_oh)
        axial_stress_oh_avg = -f_oh_sash * hoop_stress_oh_avg
    end

    if noslip
        # calc cross-sectional areas of OH and TF
        cx_surf_tf = pi * (R_tf_out^2 - R_tf_in^2)
        cx_surf_oh = pi * (R_oh_out^2 - R_oh_in^2)

        # average TF and OH axial stresses over combined TF & OH cross-sectional area, and sum
        axial_stress_tf_comb = axial_stress_tf_avg * cx_surf_tf / (cx_surf_tf + cx_surf_oh)
        axial_stress_oh_comb = axial_stress_oh_avg * cx_surf_oh / (cx_surf_tf + cx_surf_oh)
        axial_stress_comb = axial_stress_tf_comb + axial_stress_oh_comb
        axial_stress_tf_avg = axial_stress_comb
        axial_stress_oh_avg = axial_stress_comb
    end

    if plug
        r_pl = LinRange(0, R_pl, n_points)[2:end]
        displacement_pl = u_pl.(r_pl)
        ddiplacementdr_pl = dudr_pl.(r_pl)
        radial_stress_pl = sr.(r_pl, em_pl, gam_pl, displacement_pl, ddiplacementdr_pl)
        hoop_stress_pl = sh.(r_pl, em_pl, gam_pl, displacement_pl, ddiplacementdr_pl)
        axial_stress_pl_avg = 0.0
        vonmises_stress_pl = svm.(radial_stress_pl, hoop_stress_pl, axial_stress_pl_avg)
    end

    ## calculate Von Mises stress
    vonmises_stress_tf = svm.(radial_stress_tf, hoop_stress_tf, axial_stress_tf_avg)
    vonmises_stress_oh = svm.(radial_stress_oh, hoop_stress_oh, axial_stress_oh_avg)

    # apply structural composition fractions to get effective stress on structural materials
    radial_stress_tf *= 1.0 / f_struct_tf
    hoop_stress_tf *= 1.0 / f_struct_tf
    axial_stress_tf_avg *= 1.0 / f_struct_tf
    vonmises_stress_tf *= 1.0 / f_struct_tf
    radial_stress_oh *= 1.0 / f_struct_oh
    hoop_stress_oh *= 1.0 / f_struct_oh
    axial_stress_oh_avg *= 1.0 / f_struct_oh
    vonmises_stress_oh *= 1.0 / f_struct_oh
    if plug
        radial_stress_pl *= 1.0 / f_struct_pl
        hoop_stress_pl *= 1.0 / f_struct_pl
        axial_stress_pl_avg *= 1.0 / f_struct_pl
        vonmises_stress_pl *= 1.0 / f_struct_pl
    end

    smcs.bucked = Int(bucked)
    smcs.noslip = Int(noslip)
    smcs.plug = Int(plug)

    smcs.grid.r_tf = r_tf
    smcs.grid.r_oh = r_oh
    if plug
        smcs.grid.r_pl = r_pl
    end

    # keep the worse case based on the Von Mises stresses
    if ismissing(smcs.stress.vonmises, :tf) || maximum(abs, smcs.stress.vonmises.tf) < maximum(abs, vonmises_stress_tf)
        smcs.stress.vonmises.tf = vonmises_stress_tf
        smcs.stress.axial.tf = axial_stress_tf_avg
        smcs.stress.radial.tf = radial_stress_tf
        smcs.stress.hoop.tf = hoop_stress_tf
        smcs.displacement.tf = displacement_tf
    end
    if ismissing(smcs.stress.vonmises, :oh) || maximum(abs, smcs.stress.vonmises.oh) < maximum(abs, vonmises_stress_oh)
        smcs.stress.vonmises.oh = vonmises_stress_oh
        smcs.stress.axial.oh = axial_stress_oh_avg
        smcs.stress.radial.oh = radial_stress_oh
        smcs.stress.hoop.oh = hoop_stress_oh
        smcs.displacement.oh = displacement_oh
    end
    if plug && (ismissing(smcs.stress.vonmises, :pl) || maximum(abs, smcs.stress.vonmises.pl) < maximum(abs, vonmises_stress_pl))
        smcs.stress.vonmises.pl = vonmises_stress_pl
        smcs.stress.axial.pl = axial_stress_pl_avg
        smcs.stress.radial.pl = radial_stress_pl
        smcs.stress.hoop.pl = hoop_stress_pl
        smcs.displacement.pl = displacement_pl
    end

    return smcs
end
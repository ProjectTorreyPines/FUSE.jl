import Statistics

"""
    LTS_Jcrit(
        Bext,               # : strength of external magnetic field, Tesla
        strain=0.0,         # : strain on conductor due to either thermal expansion or JxB stress
        temperature=4.2,    # : temperature of conductor, Kelvin 
        )

Calculates the critical current density and fraction of critical magnetic field for internal-tin Nb3Sn LTS superconductor wire.
This uses the experimental fits from Table 6 in Lu et al., Supercond. Sci. Technol. 21 (2008) 105016.

NOTE: this returns the "engineering critical current density", which is the critical current divided by the cross-section of
the entire Nb3Sn strand. The strands considered in the reference study are approx. 50% Nb3Sn, 50% copper, and so the acutal
J_crit of the so-called "non-Cu" part of the wire (i.e. Nb3Sn only) will be approx. twice as large as the value calculated here. 

OUTPUTS
J_c : engineering critical current density, A/m^2
b   : ratio of external magnetic field at conductor to SC critical magnetic field, T/T
"""
function Nb3Sn_Jcrit(Bext, strain = 0.0, temperature = 4.2)
    #Table 6 in Lu et al., Supercond. Sci. Technol. 21 (2008) 105016
    A0 = 29330000
    Bc00 = 28.45
    Tc0 = 17.5
    c2 = -0.7388
    c3 = -0.5060
    c4 = -0.0831
    p = 0.8855
    q = 2.169
    n = 2.5
    v = 1.5
    w = 2.2
    u = 0.0
    epsilon_m = 0.0739
    epsilon = strain - epsilon_m
    Bc01 = Bc00 * (1 + c2 * epsilon^2 + c3 * epsilon^3 + c4 * epsilon^4)
    Tc = Tc0 * (Bc01 / Bc00)^(1 / w)
    t = temperature / Tc
    A = A0 * (Bc01 / Bc00)^(u / w)
    Bc11 = Bc01 * (1 - t^v)
    b = min(Bext / Bc11, 1)

    # Neutron irradiation correction based on T Baumgartner et al 2014 Supercond. Sci. Technol. 27 015005
    # if neutronFluence > 0.0
    #     p2 = 1.
    #     q2 = 2.
    #     beta_b = 1.0-exp(-(neutronFluence/1.87e22)^0.656)
    #     alpha_b = 1. - beta_b
    #     J_c = 1.0e-6*A*Tc*Tc*(1.0-t*t)^2*Bc11^(n-3)*(alpha_b*(b^(p-1)*(1-b)^q)+beta_b*(b^p2*(1-b)^q2))  # MA/m^2
    # end

    # calc critical current density
    J_c = A * Tc * Tc * (1.0 - t * t)^2 * Bc11^(n - 3) * b^(p - 1) * (1 - b)^q # A/m^2   #Equation 5 in Lu et al.

    return J_c, b
end

"""
Calculates the critical current density of YBCO high-temperature superconductor.
NOTE: the output critical current density is in terms of current per YBCO cross-sectional area. However, YBCO only comprises ~2% of the cross-sectional
area of REBCO tape. To get the critical current per cross-sectional area of REBCO tape, scale by the correct ratio of YBCO area to tape area.

    YBCO_Jcrit(
        Bext, # : external magnetic field at conductor, Tesla
        strain = 0.0, # : strain at conductor, percent
        temperature = 20.0, # : temperature of conductor, K 
        ag_c = 0) # : angle between external field and HTS tape normal vector, degrees. 0 is perp (worst-case), 90 is parallel (best-case). 

OUTPUTS
J_c : critical current density, A/m^2
b   : ratio of peak magnetic field at conductor to SC critical magnetic field, T/T
"""
function YBCO_Jcrit(Bext, strain = 0.0, temperature = 20.0, ag_c = 0)
    # Equation 2-4, Table 1 in T.S. Lee 2015 FED
    Bcrit = 68.5
    Tcrit = 87.6
    U_alp = 2.12814570695099000000E+11 # scaled to match SuperPower 2G HTS REBCO experimental results 
    # U_alp = 1.08e11    # original value from Lee 2015 FED
    V_alp = 2.0
    U_c = 0.025
    V_c = -1.2
    U_e = -0.51
    V_e = 1.09
    U_beta = 13.8
    V_beta = 0.42
    Bc2 = 0.61

    Alpha_T = U_alp * (1 - temperature / Tcrit)^V_alp
    C_T = U_c * (1 - temperature / Tcrit)^V_c
    eps0_T = U_e * (1 - temperature / Tcrit)^V_e
    beta_T = U_beta * (1 - temperature / Tcrit)^V_beta
    Bcrit_T = Bcrit * (1 - (temperature / Tcrit)^Bc2)
    eps = 1.0 - C_T * (strain - eps0_T)^2.0
    b1 = max(1 - Bext / Bcrit_T, 0.0)
    b2 = exp(-Bext * cos(ag_c * 3.1416 / 180.0) / beta_T)
    J_c = Alpha_T * eps * b1 * b2
    b = Bext / Bcrit_T

    return J_c, b
end

"""
Calculates the critical current density of ReBCO superconducting tape.

    ReBCO_Jcrit(
        Bext, # : external magnetic field at conductor, Tesla
        strain = 0.0, # : strain at conductor, percent
        temperature = 20.0, # : temperature of conductor, K 
        ag_c = 0) # : angle between external field and HTS tape normal vector, degrees. 0 is perp (worst-case), 90 is parallel (best-case). 

OUTPUTS
J_c : critical current density, A/m^2
b   : ratio of peak magnetic field at conductor to SC critical magnetic field, T/T
"""
function ReBCO_Jcrit(Bext, strain = 0.0, temperature = 20.0, ag_c = 0)
    fHTSinTape = 1.0 / 46.54 # fraction of ReBCO tape that is YBCO superconductor
    J_c, b = YBCO_Jcrit(Bext, strain, temperature, ag_c)
    return J_c * fHTSinTape, b
end

"""
    coil_technology(technology::Symbol)

Return coil parameters depending of technology [:copper, :LTS, :HTS]
"""
function coil_technology(technology::Symbol)
    coil_tech = Parameters(:coil_technology)
    if technology == :copper
        coil_tech.material = "copper"
        coil_tech.temperature = 293.0
        coil_tech.fraction_stainless = 0.0
        coil_tech.fraction_void = 0.1
    elseif technology in [:LTS, :HTS]
        if technology == :LTS
            coil_tech.temperature = 4.2
            coil_tech.material = "Nb3Sn"
        else
            coil_tech.temperature = 20.0
            coil_tech.material = "ReBCO"
        end
        coil_tech.fraction_stainless = 0.5
        coil_tech.ratio_SC_to_copper = 1.0
        coil_tech.fraction_void = 0.1
    else
        error("Supported coil tecnologies are [:copper, :LTS, :HTS")
    end
    return set_new_base!(coil_tech)
end

"""
    coil_technology(gasc::GASC, coil_type::Symbol)

Return coil parameters from GASC solution and coil type [:OH, :TF, :PF]
"""
function coil_technology(gasc::GASC, coil_type::Symbol)
    gascsol = gasc.solution
    if !(coil_type in [:OH, :TF, :PF])
        error("Supported coil type are [:OH, :TF, :PF]")
    end
    if gascsol["INPUTS"]["conductors"]["superConducting"] == "copper"
        coil_tech = coil_technology(:copper)
    else
        if gascsol["INPUTS"]["conductors"]["superConducting"] == "LTS"
            coil_tech = coil_technology(:LTS)
        elseif gascsol["INPUTS"]["conductors"]["superConducting"] == "HTS"
            coil_tech = coil_technology(:HTS)
        end
        coil_tech.thermal_strain = gascsol["INPUTS"]["conductors"]["structuralStrain$coil_type"]
        coil_tech.JxB_strain = gascsol["INPUTS"]["conductors"]["structuralStrain$coil_type"]
    end
    coil_tech.fraction_void = gascsol["INPUTS"]["conductors"]["fractionVoid$coil_type"]
    coil_tech.fraction_stainless = gascsol["INPUTS"]["conductors"]["fractionStainless$coil_type"]
    return set_new_base!(coil_tech)
end

"""
    coil_technology(machine::Symbol, coil_type::Symbol)

Return coil parameters from machine and coil type [:OH, :TF, :PF]"
"""
function coil_technology(machine::Symbol, coil_type::Symbol)
    if !(coil_type in [:OH, :TF, :PF])
        error("Supported coil type are [:OH, :TF, :PF]")
    end
    if machine == :ITER
        coil_tech = coil_technology(:LTS)
        if coil_type == :OH
            coil_tech.thermal_strain = -0.64
            coil_tech.JxB_strain = -0.05
            coil_tech.fraction_stainless = 0.46
        elseif coil_type == :TF
            coil_tech.thermal_strain = -0.69
            coil_tech.JxB_strain = -0.13
            coil_tech.fraction_stainless = 0.55
        elseif coil_type == :PF
            coil_tech.thermal_strain = -0.64
            coil_tech.JxB_strain = -0.05
            coil_tech.fraction_stainless = 0.46
        end
    else
        error("Supported coil machines are [:ITER]")
    end
    return set_new_base!(coil_tech)
end

"""
    coil_Jcrit(Bext, coil_tech::Parameters)

Returns critical current given a coil technology and external magnetic field
"""
function coil_Jcrit(Bext, coil_tech::Union{IMAS.build__pf_active__technology,IMAS.build__oh__technology,IMAS.build__tf__technology})
    if coil_tech.material == "copper"
        Jcrit = 18.5e6 # A/m^2
        fraction_cable = 1.0 - coil_tech.fraction_stainless - coil_tech.fraction_void # fraction of coil that is LTS cabling
        return Jcrit * fraction_cable # A/m^2
    else
        if coil_tech.material == "Nb3Sn"
            Jcrit_SC, Bcrit_ratio = Nb3Sn_Jcrit(Bext, coil_tech.thermal_strain + coil_tech.JxB_strain, coil_tech.temperature) # A/m^2
        elseif coil_tech.material == "ReBCO"
            Jcrit_SC, Bcrit_ratio = ReBCO_Jcrit(Bext, coil_tech.thermal_strain + coil_tech.JxB_strain, coil_tech.temperature) # A/m^2
        end
        fraction_cable = 1.0 - coil_tech.fraction_stainless - coil_tech.fraction_void # fraction of coil that is LTS cabling
        fraction_SC = fraction_cable * coil_tech.ratio_SC_to_copper / (1.0 + coil_tech.ratio_SC_to_copper) # fraction of coil that is Nb3Sn superconductor
        return Jcrit_SC * fraction_SC # A/m^2
    end
end


"""
    solve_1D_solid_mechanics(
        R0,                 # : (float) major radius at center of TF bore, meters
        B0,                 # : (float) toroidal field at R0, Tesla
        R_tf_in,            # : (float) major radius of inboard edge of TF coil core legs, meters
        R_tf_out,           # : (float) major radius of outboard edge of TF coil core legs, meters
        Bz_cs,              # : (float) axial field in solenoid bore, Tesla
        R_cs_in,            # : (float) major radius of inboard edge of CS coil, meters
        R_cs_out,           # : (float) major radius of outboard edge of CS coil, meters
        sz_tf_avg = nothing,   # : (float) average axial stress in TF coil core legs, Pa (if nothing, use constant fraction of hoop stress)
        sz_cs_avg = nothing,   # : (float) average axial stress in CS coil, Pa (if nothing, use constant fraction of hoop stress)
        TFCSbucked = false, # : (bool), flag for bucked boundary conditions between TF and CS (and center plug, if present)
        noslip = false,     # : (bool), flag for no slip conditions between TF and CS (and center plug, if present)
        doplug = false,     # : (bool), flag for center plug
        f_struct_tf = 1.0,   # : (float), fraction of TF coil that is structural material
        f_struct_cs = 1.0,   # : (float), fraction of CS coil that is structural material
        f_struct_pl = 1.0,   # : (float), fraction of plug that is structural material
        em_tf = 193103448275.0,# : (float), modulus of elasticity for TF coil, Pa (default is stainless steel)
        gam_tf = 0.33,        # : (float), Poisson"s ratio for TF coil, (default is stainless steel)
        em_cs = 193103448275.0,# : (float), modulus of elasticity for CS coil, Pa (default is stainless steel)
        gam_cs = 0.33,        # : (float), Poisson"s ratio for CS coil, (default is stainless steel)
        em_pl = 193103448275.0,# : (float), modulus of elasticity for center plug, Pa (default is stainless steel)
        gam_pl = 0.33,        # : (float), Poisson"s ratio for center plug, (default is stainless steel)
        f_tf_sash = 0.873,  # : (float), conversion factor from hoop stress to axial stress for TF coil (nominally 0.873)
        f_cs_sash = 0.37337,  # : (float), conversion factor from hoop stress to axial stress for CS coil (nominally 0.37337)
        verbose = false,      # : (bool), flag for verbose output to terminal
        return_profs = false, # : (bool), flag to include profiles of stress in output
)

Uses Leuer 1D solid mechanics equations to solve for radial and hoop stresses in TF coil, CS coil, and center plug.
Based on derivations in Engineering Physics Note "EPNjal17dec17_gasc_pt5_tf_cs_plug_buck" by Jim Leuer (Dec. 17, 2017)

Returns radial, hoop, axial, and Von Mises stresses (average and peak) for TF, CS, and plug (Pascals)
(optional) radial profiles of radial, hoop, axial, and Von Mises stresses for TF, CS, and plug (Pascals)

The tokamak radial buid is :

|| plug or void (0 < r < R1) || coil 1 (R1 < r < R2) || coil 2 (R3 < r < R4) || ----> plasma center (r = R0)

NOTE: effect of structural composition fractions < 1 is only applied to Von Mises stress (not radial, hoop, or axial)!!
"""
function solve_1D_solid_mechanics(
    R0,                    # : (float) major radius at center of TF bore, meters
    B0,                    # : (float) toroidal field at R0, Tesla
    R_tf_in,               # : (float) major radius of inboard edge of TF coil core legs, meters
    R_tf_out,              # : (float) major radius of outboard edge of TF coil core legs, meters
    Bz_cs,                 # : (float) axial field in solenoid bore, Tesla
    R_cs_in,               # : (float) major radius of inboard edge of CS coil, meters
    R_cs_out;              # : (float) major radius of outboard edge of CS coil, meters
    sz_tf_avg = nothing,   # : (float) average axial stress in TF coil core legs, Pa (if nothing, use constant fraction of hoop stress)
    sz_cs_avg = nothing,   # : (float) average axial stress in CS coil, Pa (if nothing, use constant fraction of hoop stress)
    TFCSbucked = false,    # : (bool), flag for bucked boundary conditions between TF and CS (and center plug, if present)
    noslip = false,        # : (bool), flag for no slip conditions between TF and CS (and center plug, if present)
    doplug = false,        # : (bool), flag for center plug
    f_struct_tf = 1.0,     # : (float), fraction of TF coil that is structural material
    f_struct_cs = 1.0,     # : (float), fraction of CS coil that is structural material
    f_struct_pl = 1.0,     # : (float), fraction of plug that is structural material
    em_tf = 193103448275.0,# : (float), modulus of elasticity for TF coil, Pa (default is stainless steel)
    gam_tf = 0.33,         # : (float), Poisson"s ratio for TF coil, (default is stainless steel)
    em_cs = 193103448275.0,# : (float), modulus of elasticity for CS coil, Pa (default is stainless steel)
    gam_cs = 0.33,         # : (float), Poisson"s ratio for CS coil, (default is stainless steel)
    em_pl = 193103448275.0,# : (float), modulus of elasticity for center plug, Pa (default is stainless steel)
    gam_pl = 0.33,         # : (float), Poisson"s ratio for center plug, (default is stainless steel)
    f_tf_sash = 0.873,     # : (float), conversion factor from hoop stress to axial stress for TF coil (nominally 0.873)
    f_cs_sash = 0.37337,   # : (float), conversion factor from hoop stress to axial stress for CS coil (nominally 0.37337)
    verbose = false        # : (bool), flag for verbose output to terminal
)

    if verbose
        println("solve_1D_solid_mechanics:")
        println("- R0 = ", R0)
        println("- B0 = ", B0)
        println("- R_tf_in = ", R_tf_in)
        println("- R_tf_out = ", R_tf_out)
        println("- Bz_cs = ", Bz_cs)
        println("- R_cs_in = ", R_cs_in)
        println("- R_cs_out = ", R_cs_out)
        println("- sz_tf_avg = ", sz_tf_avg)
        println("- sz_cs_avg = ", sz_cs_avg)
        println("- TFCSbucked = ", TFCSbucked)
        println("- noslip = ", noslip)
        println("- doplug = ", doplug)
        println("- f_struct_tf", f_struct_tf)
        println("- f_struct_cs", f_struct_cs)
        println("- f_struct_pl", f_struct_pl)
    end

    # define structural constants
    mu0 = 4e-7 * pi
    embar_tf = em_tf / (1 - gam_tf^2)
    embar_cs = em_cs / (1 - gam_cs^2)
    embar_pl = em_pl / (1 - gam_pl^2)

    # define forcing constraints on TF and CS coils
    C_tf = 1.0 / embar_tf * 2 * (B0 * R0)^2 / (mu0 * (R_tf_out^2 - R_tf_in^2)^2)
    C_cs = -1.0 / embar_cs * Bz_cs^2 / (mu0 * (R_cs_out - R_cs_in)^2)

    # calc centerlines, check radial build inputs for consistency 
    cl_tf = 0.5 * (R_tf_out + R_tf_in)
    cl_cs = 0.5 * (R_cs_out + R_cs_in)

    # determine ordering of radial build
    if cl_tf < cl_cs
        order = "TF-CS"
    else
        order = "CS-TF"
    end
    if verbose
        println("* order = ", order)
    end

    # check radial build inputs for consistency 
    if R_tf_out < R_tf_in
        error("R_tf_out ({:.3f}) < R_tf_in ({:.3f})".format(R_tf_out, R_tf_in))
    end
    if R_cs_out < R_cs_in
        error("R_cs_out ({:.3f}) < R_cs_in ({:.3f})".format(R_cs_out, R_cs_in))
    end
    # order == "TF-CS" and R_cs_in < R_tf_out
    if (order == "TF-CS") && ((R_tf_out - R_cs_in * (1 + 1e-3)) > 0)
        println("* TF and CS are overlapping: R_cs_in ($R_cs_in) < R_tf_out ($R_tf_out)")
    end
    # order == "CS-TF" and R_tf_in < R_cs_out
    if (order == "CS-TF") && ((R_cs_out - R_tf_in * (1 + 1e-3)) > 0)
        println("* TF and CS are overlapping: R_tf_in ($R_tf_in) < R_cs_out ($R_cs_out)")
    end
    # TFCSbucked and order == "CS-TF" and R_cs_out != R_tf_in
    if TFCSbucked && (order == "CS-TF") && (abs.(R_cs_out - R_tf_in) > R_tf_in * 1e-3)
        println("* CS is bucked against TF, but R_cs_out ($R_cs_out) != R_tf_in ($R_tf_in)")
    end

    # define radial functions in u/r and du/dr for TF and CS
    # u_tf/r   = C_tf * f_tf(r) + A_tf + B_tf/r^2
    # du_tf/dr = C_tf * g_tf(r) + A_tf - B_tf/r^2
    # u_cs/r   = C_cs * f_cs(r) + A_cs + B_cs/r^2
    # du_cs/dr = C_cs * g_cs(r) + A_cs - B_cs/r^2
    function f_tf(r)
        logr = log(r)
        try
            logr[!isfinite(logr)] = -100.0 # avoid nan output for non-positive r
        catch
            if !isfinite(logr)
                logr = -100.0
            end
            return 1.0 / 8.0 * r^2 - R_tf_in^2 / 2.0 * (logr - 0.5)
        end
    end

    function g_tf(r)
        logr = log(r)
        try
            logr[!isfinite(logr)] = -100.0 # avoid nan output for non-positive r
        catch
            if !isfinite(logr)
                logr = -100.0
            end
            return 3.0 / 8.0 * r^2 - R_tf_in^2 / 2.0 * (logr + 0.5)
        end
    end

    function f_cs(r)
        return 1.0 / 3.0 * R_cs_out * r - 1.0 / 8.0 * r^2
    end
    function g_cs(r)
        return 2.0 / 3.0 * R_cs_out * r - 3.0 / 8.0 * r^2
    end

    # define linear system of equations based on boundary conditions
    # M * X = Y

    ## free standing CS and TF coils
    # radial stress = 0 at ALL coil edges
    if !TFCSbucked
        doplug = false
        if verbose
            println("* Free standing coils")
        end

        M = [
            [1 + gam_tf, (gam_tf - 1) / R_tf_in^2, 0.0, 0.0];;
            [1 + gam_tf, (gam_tf - 1) / R_tf_out^2, 0.0, 0.0];;
            [0.0, 0.0, 1 + gam_cs, (gam_cs - 1) / R_cs_in^2];;
            [0.0, 0.0, 1 + gam_cs, (gam_cs - 1) / R_cs_out^2]
        ]
        Y = [-C_tf * (g_tf(R_tf_in) + gam_tf * f_tf(R_tf_in)),
            -C_tf * (g_tf(R_tf_out) + gam_tf * f_tf(R_tf_out)),
            -C_cs * (g_cs(R_cs_in) + gam_cs * f_cs(R_cs_in)),
            -C_cs * (g_cs(R_cs_out) + gam_cs * f_cs(R_cs_out)),
        ]
        A_tf, B_tf, A_cs, B_cs = M \ Y
        A_pl = 0.0

        ## bucked CS and TF only (no plug)
        # order must be "CS-TF"
        # radial stress = 0 at R_cs_in and R_tf_out
        # radial stresses and displacementes are equal at interface R_int ( == R_cs_out = R_tf_in )
    elseif !doplug && (order == "CS-TF")
        if verbose
            println("* CS bucked against TF only")
        end
        R_int = R_cs_out
        M = [
            [1 + gam_tf, (gam_tf - 1) / R_tf_out^2, 0.0, 0.0];;
            [0.0, 0.0, 1 + gam_cs, (gam_cs - 1) / R_cs_in^2];;
            [1.0, 1.0 / R_int^2, -1.0, -1.0 / R_int^2];;
            [embar_tf * (1 + gam_tf), embar_tf * (gam_tf - 1) / R_int^2, -embar_cs * (1 + gam_cs), -embar_cs * (gam_cs - 1) / R_int^2]
        ]
        Y = [-C_tf * (g_tf(R_tf_out) + gam_tf * f_tf(R_tf_out)),
            -C_cs * (g_cs(R_cs_in) + gam_cs * f_cs(R_cs_in)),
            C_cs * f_cs(R_int) - C_tf * f_tf(R_int),
            embar_cs * C_cs * (g_cs(R_int) + gam_cs * f_cs(R_int)) - embar_tf * C_tf * (g_tf(R_int) + gam_tf * f_tf(R_int)),
        ]
        A_tf, B_tf, A_cs, B_cs = M \ Y
        A_pl = 0.0

        ## bucked plug, CS, and TF
        # order is "CS-TF"
        # radial stress = 0 at R_tf_out
        # radial stresses and displacements are equal at plug-CS interface R_pl ( == R_cs_in) 
        #  and CS-TF interface R_int ( == R_cs_out = R_tf_in )
    elseif doplug && (order == "CS-TF")
        if verbose
            println("* CS bucked against TF and plug")
        end
        R_int = R_cs_out
        R_pl = R_cs_in
        M = [
            [1 + gam_tf, (gam_tf - 1) / R_tf_out^2, 0.0, 0.0, 0.0];;
            [1.0, 1.0 / R_int^2, -1.0, -1.0 / R_int^2, 0.0];;
            [embar_tf * (1 + gam_tf), embar_tf * (gam_tf - 1) / R_int^2, -embar_cs * (1 + gam_cs), -embar_cs * (gam_cs - 1) / R_int^2, 0.0];;
            [0.0, 0.0, 1.0, 1.0 / R_pl^2, -1.0];;
            [0.0, 0.0, embar_cs * (1 + gam_cs), embar_cs * (gam_cs - 1) / R_pl^2, -embar_pl * (1 + gam_pl)]
        ]
        Y = [-C_tf * (g_tf(R_tf_out) + gam_tf * f_tf(R_tf_out)),
            C_cs * f_cs(R_int) - C_tf * f_tf(R_int),
            embar_cs * C_cs * (g_cs(R_int) + gam_cs * f_cs(R_int)) - embar_tf * C_tf * (g_tf(R_int) + gam_tf * f_tf(R_int)),
            -C_cs * f_cs(R_pl),
            -embar_cs * C_cs * (g_cs(R_pl) + gam_cs * f_cs(R_pl)),
        ]
        A_tf, B_tf, A_cs, B_cs, A_pl = M \ Y

    ## bucked TF aginst plug, CS is free standing
    # order is "TF-CS"
    # radial stress = 0 at R_tf_out
    # radial stress = 0 at R_cs_in and R_cs_out
    # radial stresses and displacements are equal at plug-TF interface R_pl ( == R_tf_in)
    elseif doplug && (order == "TF-CS")
        if verbose
            println("* TF bucked against plug")
        end
        R_pl = R_tf_in
        M = [
            [1 + gam_tf, (gam_tf - 1) / R_tf_out^2, 0.0, 0.0, 0.0];;
            [0.0, 0.0, 1 + gam_cs, (gam_cs - 1) / R_cs_in^2, 0.0];;
            [0.0, 0.0, 1 + gam_cs, (gam_cs - 1) / R_cs_out^2, 0.0];;
            [1.0, 1.0 / R_pl^2, 0.0, 0.0, -1.0];;
            [embar_tf * (1 + gam_tf), embar_tf * (gam_tf - 1) / R_pl^2, 0.0, 0.0, -embar_pl * (1 + gam_pl)]
        ]
        Y = [-C_tf * (g_tf(R_tf_out) + gam_tf * f_tf(R_tf_out)),
            -C_cs * (g_cs(R_cs_in) + gam_cs * f_cs(R_cs_in)),
            -C_cs * (g_cs(R_cs_out) + gam_cs * f_cs(R_cs_out)),
            -C_tf * f_tf(R_pl),
            -embar_tf * C_tf * (g_tf(R_pl) + gam_tf * f_tf(R_pl)),
        ]
        A_tf, B_tf, A_cs, B_cs, A_pl = M \ Y

    else
        error("radial build inputs do not match known boundary conditions!!")
    end

    ## use integration constants to define functions for displacement and radial, hoop, and Vom Mises stresses

    function u_tf(r)
        return @. r * (C_tf * f_tf(r) + A_tf + B_tf / r^2)
    end

    function dudr_tf(r)
        return @. C_tf * g_tf(r) + A_tf - B_tf / r^2
    end

    function u_cs(r)
        return @. r * (C_cs * f_cs(r) + A_cs + B_cs / r^2)
    end

    function dudr_cs(r)
        return @. C_cs * g_cs(r) + A_cs - B_cs / r^2
    end

    function u_pl(r)
        return @. r * A_pl
    end

    function dudr_pl(r)
        return @. A_pl
    end

    function sr(r, em, gam, u, dudr)
        return @. em / (1 - gam^2) * (dudr + u / r * gam)
    end

    function sh(r, em, gam, u, dudr)
        return @. em / (1 - gam^2) * (u / r + dudr * gam)
    end

    function svm(sr, sh, sa)
        return @. sqrt(((sh - sa)^2 + (sa - sr)^2 + (sr - sh)^2) / 2.0)
    end

    ## calc radial and hoop stresses (average and peak)
    # also estimate axial stresses if not given & modify due to noslip condition
    # default axial scalings taken from Bending free formula & GASC
    nr = 51
    r_cs = LinRange(R_cs_in, R_cs_out, nr)
    r_tf = LinRange(R_tf_in, R_tf_out, nr)
    u_cs_arr = u_cs(r_cs)
    u_tf_arr = u_tf(r_tf)
    sr_cs_arr = sr(r_cs, em_cs, gam_cs, u_cs(r_cs), dudr_cs(r_cs))
    sr_tf_arr = sr(r_tf, em_tf, gam_tf, u_tf(r_tf), dudr_tf(r_tf))
    sh_cs_arr = sh(r_cs, em_cs, gam_cs, u_cs(r_cs), dudr_cs(r_cs))
    sh_tf_arr = sh(r_tf, em_tf, gam_tf, u_tf(r_tf), dudr_tf(r_tf))

    sr_cs_avg = Statistics.mean(sr_cs_arr)
    sr_tf_avg = Statistics.mean(sr_tf_arr)
    sh_cs_avg = Statistics.mean(sh_cs_arr)
    sh_tf_avg = Statistics.mean(sh_tf_arr)

    sr_cs_peak = maximum(abs.(sr_cs_arr))
    sr_tf_peak = maximum(abs.(sr_tf_arr))
    sh_cs_peak = maximum(abs.(sh_cs_arr))
    sh_tf_peak = maximum(abs.(sh_tf_arr))

    if sz_tf_avg === nothing
        sz_tf_avg0 = -f_tf_sash * sh_tf_avg
    else
        sz_tf_avg0 = sz_tf_avg
    end
    if sz_cs_avg === nothing
        sz_cs_avg0 = -f_cs_sash * sh_cs_avg
    else
        sz_cs_avg0 = sz_cs_avg
    end

    if noslip
        # calc cross-sectional areas of CS and TF
        cx_surf_tf = pi * (R_tf_out^2 - R_tf_in^2)
        cx_surf_cs = pi * (R_cs_out^2 - R_cs_in^2)

        # average TF and CS axial stresses over combined TF & CS cross-sectional area, and sum
        sz_tf_comb = sz_tf_avg0 * cx_surf_tf / (cx_surf_tf + cx_surf_cs)
        sz_cs_comb = sz_cs_avg0 * cx_surf_cs / (cx_surf_tf + cx_surf_cs)
        sz_comb = sz_tf_comb + sz_cs_comb
        sz_tf_avg0 = sz_comb
        sz_cs_avg0 = sz_comb
    end

    if doplug
        r_pl = LinRange(0, R_pl, nr)[2:end]
        u_pl_arr = u_pl(r_pl)
        sr_pl_arr = sr(r_pl, em_pl, gam_pl, u_pl(r_pl), dudr_pl(r_pl))
        sh_pl_arr = sh(r_pl, em_pl, gam_pl, u_pl(r_pl), dudr_pl(r_pl))
        sr_pl_avg = Statistics.mean(sr_pl_arr)
        sh_pl_avg = Statistics.mean(sh_pl_arr)
        sr_pl_peak = maximum(abs.(sr_pl_arr))
        sh_pl_peak = maximum(abs.(sh_pl_arr))
        sz_pl_avg = 0.0
    else
        r_pl = [0.0]
        u_pl_arr = [0.0]
        sr_pl_arr = [0.0]
        sh_pl_arr = [0.0]
        sr_pl_avg = 0.0
        sh_pl_avg = 0.0
        sr_pl_peak = 0.0
        sh_pl_peak = 0.0
        sz_pl_avg = 0.0
    end

    ## calculate Von Mises stress (average and peak)
    svm_tf_arr = svm(sr_tf_arr, sh_tf_arr, sz_tf_avg0)
    svm_cs_arr = svm(sr_cs_arr, sh_cs_arr, sz_cs_avg0)

    svm_tf_avg = Statistics.mean(svm_tf_arr)
    svm_cs_avg = Statistics.mean(svm_cs_arr)

    svm_tf_peak = maximum(abs.(svm_tf_arr))
    svm_cs_peak = maximum(abs.(svm_cs_arr))

    if doplug
        svm_pl_arr = svm(sr_pl_arr, sh_pl_arr, sz_pl_avg)
        svm_pl_avg = Statistics.mean(svm_cs_arr)
        svm_pl_peak = maximum(abs.(svm_cs_arr))
    else
        svm_pl_arr = [0.0]
        svm_pl_avg = 0.0
        svm_pl_peak = 0.0
    end

    # apply structural composition fractions to get effective Von Mises stress on structural materials
    sr_tf_arr *= 1.0 / f_struct_tf
    sr_tf_avg *= 1.0 / f_struct_tf
    sr_tf_peak *= 1.0 / f_struct_tf
    sh_tf_arr *= 1.0 / f_struct_tf
    sh_tf_avg *= 1.0 / f_struct_tf
    sh_tf_peak *= 1.0 / f_struct_tf
    sz_tf_avg0 *= 1.0 / f_struct_tf
    svm_tf_arr *= 1.0 / f_struct_tf
    svm_tf_avg *= 1.0 / f_struct_tf
    svm_tf_peak *= 1.0 / f_struct_tf

    sr_cs_arr *= 1.0 / f_struct_cs
    sr_cs_avg *= 1.0 / f_struct_cs
    sr_cs_peak *= 1.0 / f_struct_cs
    sh_cs_arr *= 1.0 / f_struct_cs
    sh_cs_avg *= 1.0 / f_struct_cs
    sh_cs_peak *= 1.0 / f_struct_cs
    sz_cs_avg0 *= 1.0 / f_struct_cs
    svm_cs_arr *= 1.0 / f_struct_cs
    svm_cs_avg *= 1.0 / f_struct_cs
    svm_cs_peak *= 1.0 / f_struct_cs

    sr_pl_arr *= 1.0 / f_struct_pl
    sr_pl_avg *= 1.0 / f_struct_pl
    sr_pl_peak *= 1.0 / f_struct_pl
    sh_pl_arr *= 1.0 / f_struct_pl
    sh_pl_avg *= 1.0 / f_struct_pl
    sh_pl_peak *= 1.0 / f_struct_pl
    sz_pl_avg *= 1.0 / f_struct_pl
    svm_pl_arr *= 1.0 / f_struct_pl
    svm_pl_avg *= 1.0 / f_struct_pl
    svm_pl_peak *= 1.0 / f_struct_pl

    return Dict(
        "svm_tf_peak" => svm_tf_peak, "svm_tf_avg" => svm_tf_avg,
        "svm_cs_peak" => svm_cs_peak, "svm_cs_avg" => svm_cs_avg,
        "svm_pl_peak" => svm_pl_peak, "svm_pl_avg" => svm_pl_avg, "sr_tf_peak" => sr_tf_peak, "sr_tf_avg" => sr_tf_avg,
        "sr_cs_peak" => sr_cs_peak, "sr_cs_avg" => sr_cs_avg,
        "sr_pl_peak" => sr_pl_peak, "sr_pl_avg" => sr_pl_avg, "sh_tf_peak" => sh_tf_peak, "sh_tf_avg" => sh_tf_avg,
        "sh_cs_peak" => sh_cs_peak, "sh_cs_avg" => sh_cs_avg,
        "sh_pl_peak" => sh_pl_peak, "sh_pl_avg" => sh_pl_avg, "sz_tf_avg" => sz_tf_avg0,
        "sz_cs_avg" => sz_cs_avg0,
        "sz_pl_acg" => sz_pl_avg, "r_tf" => r_tf, "r_cs" => r_cs, "r_pl" => r_pl,
        "u_tf_arr" => u_tf_arr, "u_cs_arr" => u_cs_arr, "u_pl_arr" => u_pl_arr,
        "svm_tf_arr" => svm_tf_arr, "svm_cs_arr" => svm_cs_arr, "svm_pl_arr" => svm_pl_arr,
        "sr_tf_arr" => sr_tf_arr, "sr_cs_arr" => sr_cs_arr, "sr_pl_arr" => sr_pl_arr,
        "sh_tf_arr" => sh_tf_arr, "sh_cs_arr" => sh_cs_arr, "sh_pl_arr" => sh_pl_arr,
    )
end
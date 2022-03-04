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
    elseif technology in [:LTS, :HTS]
        if technology == :LTS
            coil_tech.temperature = 4.2
            coil_tech.material = "Nb3Sn"
        else
            coil_tech.temperature = 20.0
            coil_tech.material = "ReBCO"
        end
        coil_tech.ratio_SC_to_copper = 1.0
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
        Jcrit = 100000.0 # A/m^2
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
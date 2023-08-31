const supported_coils_techs = [:copper, :LTS, :ITER, :HTS]

"""
Material properties
"""
Base.@kwdef struct MaterialProperties
    yield_strength::Float64 = NaN
    young_modulus::Float64 = NaN
    poisson_ratio::Float64 = NaN
end

const stainless_steel = MaterialProperties(
    yield_strength=800E6, # Pa
    young_modulus=193.103448275E9, # Pa
    poisson_ratio=0.33,
)

const pure_copper = MaterialProperties(
    yield_strength=70E6, # Pa
    young_modulus=110E9, # Pa
    poisson_ratio=0.34,
)

"""
    coil_technology(technology::Symbol, coil_type::Symbol)

Return coil parameters from technology $(supported_coils_techs) and coil type [:oh, :tf, :pf_active]"
"""
function coil_technology(coil_tech::Union{IMAS.build__pf_active__technology,IMAS.build__oh__technology,IMAS.build__tf__technology}, technology::Symbol, coil_type::Symbol)
    if coil_type ∉ (:oh, :tf, :pf_active)
        error("Supported coil type are [:oh, :tf, :pf_active]")
    end
    if technology ∉ supported_coils_techs
        error("Supported coil technology are $(repr(supported_coils_techs))")
    end

    if technology == :copper
        coil_tech.material = "Copper"
        coil_tech.temperature = 293.0
        coil_tech.fraction_steel = 0.0
        coil_tech.fraction_void = 0.1

    elseif technology ∈ (:LTS, :ITER, :HTS)
        if technology ∈ (:LTS, :ITER)
            coil_tech.temperature = 4.2
            coil_tech.material = "Nb3Sn"
        else
            coil_tech.temperature = 4.2
            coil_tech.material = "ReBCO"
        end
        coil_tech.fraction_steel = 0.5
        coil_tech.ratio_SC_to_copper = 1.0
        coil_tech.fraction_void = 0.1
    end

    if technology == :ITER
        if coil_type == :oh
            coil_tech.thermal_strain = -0.64
            coil_tech.JxB_strain = -0.05
            coil_tech.fraction_steel = 0.46
        elseif coil_type == :tf
            coil_tech.thermal_strain = -0.69
            coil_tech.JxB_strain = -0.13
            coil_tech.fraction_steel = 0.55
        elseif coil_type == :pf_active
            coil_tech.thermal_strain = -0.64
            coil_tech.JxB_strain = -0.05
            coil_tech.fraction_steel = 0.46
        end
    end

    coil_tech.thermal_strain = 0.0
    coil_tech.JxB_strain = 0.0

    return coil_tech
end

"""
    LTS_Jcrit(
        Bext::Real,               # : strength of external magnetic field, Tesla
        strain::Real=0.0,         # : strain on conductor due to either thermal expansion or JxB stress
        temperature::Real=4.2,    # : temperature of conductor, Kelvin 
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
function Nb3Sn_Jcrit(Bext::Real, strain::Real=0.0, temperature::Real=4.2)
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
        Bext::Real, # : external magnetic field at conductor, Tesla
        strain::Real = 0.0, # : strain at conductor, percent
        temperature::Real = 20.0, # : temperature of conductor, K 
        ag_c::Real = 0.0) # : angle between external field and HTS tape normal vector, degrees. 0 is perp (worst-case), 90 is parallel (best-case). 

OUTPUTS
J_c : critical current density, A/m^2
b   : ratio of peak magnetic field at conductor to SC critical magnetic field, T/T
"""
function YBCO_Jcrit(Bext::Real, strain::Real=0.0, temperature::Real=20.0, ag_c::Real=0.0)
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
        Bext::Real, # : external magnetic field at conductor, Tesla
        strain::Real = 0.0, # : strain at conductor, percent
        temperature::Real = 20.0, # : temperature of conductor, K 
        ag_c::Real = 0.0) # : angle between external field and HTS tape normal vector, degrees. 0 is perp (worst-case), 90 is parallel (best-case). 

OUTPUTS
J_c : critical current density, A/m^2
b   : ratio of peak magnetic field at conductor to SC critical magnetic field, T/T
"""
function ReBCO_Jcrit(Bext::Real, strain::Real=0.0, temperature::Real=20.0, ag_c::Real=0.0)
    fHTSinTape = 1.0 / 46.54 # fraction of ReBCO tape that is YBCO superconductor
    J_c, b = YBCO_Jcrit(Bext, strain, temperature, ag_c)
    return J_c * fHTSinTape, b
end

"""
    coil_J_B_crit(Bext, coil_tech::Union{IMAS.build__pf_active__technology,IMAS.build__oh__technology,IMAS.build__tf__technology})

Returns critical current density and magnetic field given an external magnetic field and coil technology
"""
function coil_J_B_crit(Bext, coil_tech::Union{IMAS.build__pf_active__technology,IMAS.build__oh__technology,IMAS.build__tf__technology})
    fraction_conductor = 1.0 - coil_tech.fraction_steel - coil_tech.fraction_void # fraction of coil that is a conductor
    @assert fraction_conductor > 0.0 "coil_J_B_crit: coil technology has no room for conductor"
    if coil_tech.material == "Copper"
        Jcrit = 18.5e6 # A/m^2
        return Jcrit * fraction_conductor, Inf # A/m^2
    else
        if coil_tech.material == "Nb3Sn"
            Jcrit_SC, Bext_Bcrit_ratio = Nb3Sn_Jcrit(Bext, coil_tech.thermal_strain + coil_tech.JxB_strain, coil_tech.temperature) # A/m^2
        elseif coil_tech.material == "ReBCO"
            Jcrit_SC, Bext_Bcrit_ratio = ReBCO_Jcrit(Bext, coil_tech.thermal_strain + coil_tech.JxB_strain, coil_tech.temperature) # A/m^2
        end
        fraction_SC = fraction_conductor * coil_tech.ratio_SC_to_copper / (1.0 + coil_tech.ratio_SC_to_copper) # fraction of coil that is Nb3Sn superconductor
        Jcrit = Jcrit_SC * fraction_SC # A/m^2
        return Jcrit, Bext / Bext_Bcrit_ratio
    end
end

function GAMBL_blanket(bm::IMAS.blanket__module)
    layers = resize!(bm.layer, 3)

    n = 1
    layers[n].name = "First wall"
    layers[n].material = "Tungsten"
    layers[n].thickness = 0.02

    n = n + 1
    layers[n].name = "Breeder"
    layers[n].material = "lithium-lead"
    layers[n].thickness = 0.5

    n = n + 1
    layers[n].name = "Shield"
    layers[n].material = "Tungsten"
    layers[n].thickness = 0.05

    return bm
end

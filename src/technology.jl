const supported_coils_techs = [:copper, :Nb3Sn, :NbTi, :ITER, :KDEMO, :HTS]

"""
Material properties
"""
Base.@kwdef struct MaterialProperties
    yield_strength::Float64 = NaN
    young_modulus::Float64 = NaN
    poisson_ratio::Float64 = NaN
end

const stainless_steel = MaterialProperties(;
    yield_strength=800E6, # Pa
    young_modulus=193.103448275E9, # Pa
    poisson_ratio=0.33
)

const pure_copper = MaterialProperties(;
    yield_strength=70E6, # Pa
    young_modulus=110E9, # Pa
    poisson_ratio=0.34
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

    elseif technology ∈ (:Nb3Sn, :NbTi, :ITER, :KDEMO, :HTS)
        if technology ∈ (:Nb3Sn, :ITER)
            coil_tech.temperature = 4.2
            coil_tech.material = "Nb3Sn"
            coil_tech.fraction_void = 0.1
        elseif technology == :NbTi 
            coil_tech.temperature =  4.2
            coil_tech.material = "NbTi"
            coil_tech.fraction_void = 0.2 # from Supercond. Sci. Technol. 36 (2023) 075009
        elseif technology == :KDEMO
            coil_tech.temperature = 4.2 
            coil_tech.material = "KDEMO Nb3Sn"
            if coil_type == :tf
                coil_tech.fraction_void = 0.26 # from NF 55 (2015) 053027, Table 2
            end
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

mutable struct LTS_scaling
    A0::Real
    Bc00::Real
    epsilon_m::Real
    Tc0::Real
    c2::Real
    c3::Real
    c4::Real
    p::Real
    q::Real
    n::Real
    u::Real
    v::Real
    w::Real
end

"""
LTS_Jcrit(
    lts::LTS_scaling          # : object containing scaling parameters, specific to conductor type
    Bext::Real,               # : strength of external magnetic field, Tesla
    strain::Real=0.0,         # : strain on conductor due to either thermal expansion or JxB stress
    temperature::Real=4.2,    # : temperature of conductor, Kelvin 
)

Calculates the critical current density and fraction of critical magnetic field for either NbTi or Nb3Sn wire, depending on the parameters 
passed in the lts object. 

NOTE: this returns the "engineering critical current density", which is the critical current divided by the cross-section of
the entire Nb3Sn strand. The strands considered in the reference study are approx. 50% Nb3Sn, 50% copper, and so the acutal
J_crit of the so-called "non-Cu" part of the wire (i.e. Nb3Sn only) will be approx. twice as large as the value calculated here.

OUTPUTS
J_c : engineering critical current density, A/m^2
b   : ratio of external magnetic field at conductor to SC critical magnetic field, T/T
"""

function LTS_Jcrit(lts::LTS_scaling, Bext::Real, strain::Real=0.0, temperature::Real=4.2)
    epsilon = strain - lts.epsilon_m
    Bc01 = lts.Bc00 * (1 + lts.c2 * epsilon^2 + lts.c3 * epsilon^3 + lts.c4 * epsilon^4)
    Tc = lts.Tc0 * (Bc01 / lts.Bc00)^(1 / lts.w)
    t = temperature / Tc
    A = lts.A0 * (Bc01 / lts.Bc00)^(lts.u / lts.w)
    Bc11 = Bc01 * (1 - t^lts.v)
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
    J_c = A * Tc * Tc * (1.0 - t * t)^2 * Bc11^(lts.n - 3) * b^(lts.p - 1) * (1 - b)^lts.q # A/m^2   #Equation 5 in Lu et al.

    return J_c, b
end

"""

KDEMO_Nb3Sn_Jcrit(
    Bext::Real,         # External magnetic field at conductor, Tesla 
    strain::Real,       # Strain at conductor, percent 
    temperature::Real   # Temperature at conductor, K
)

Calculates critical current density for Nb3Sn conductor according to the ITER 2008 parametrization. The nine free parameters 
(p, q, C, Ca1, Ca2, epsilon_0a, epsilon_m, Bc20m, Tc0m) have to be determined experimentally and vary between conductors produced 
by different manufacturers. 

The specific values of the parameters here were determined experimentally for the K-DEMO Nb3Sn conductors. 

"""

function KDEMO_Nb3Sn_Jcrit(Bext::Real, strain::Real, temperature::Real=4.2)
    strain = strain ./ 1e2 # convert from percent to (mm/mm)

    # Nine free parameters determined experimentally for KDEMO superconductor 
    p = 0.5
    q = 2
    C = 67981.2 * 1e6
    Ca1 = 56.7401
    Ca2 = 22.7938
    epsilon_0a = 0.00420247
    epsilon_m = 0.00363043
    Bc20m = 28.9919
    Tc0m = 16.8

    kdemo_diameter = 0.82 # mm 
    kdemo_copper_frac = 0.5 # strands are 50% copper, 50% Nb3Sn

    epsilon = strain - epsilon_m 
    epsilon_sh = (Ca2*epsilon_0a)/(sqrt((Ca1)^2 - (Ca2)^2))

    s1 = Ca1 * (sqrt(epsilon_sh^2 + epsilon_0a^2) - sqrt((epsilon - epsilon_sh)^2 + epsilon_0a^2))
    s2 = Ca2 * epsilon
    s3 = 1 - Ca1*epsilon_0a
    
    s = 1 + ((s1 - s2)/s3)

    Tc_star = Tc0m * s^(1/3)
    t = temperature/Tc_star 

    Bc_star = Bc20m * s * (1 - t^1.52)
    b = Bext/Bc_star

    Ic = (C/Bext) * s * (1.0 - t^1.52) * (1.0 - t^2) * b^p * ((1.0 - b)^q)

    Jc = Ic / (pi*(0.5 * kdemo_diameter)^2 * kdemo_copper_frac)

    return Jc, b
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
            # Jcrit_SC, Bext_Bcrit_ratio = Nb3Sn_Jcrit(Bext, coil_tech.thermal_strain + coil_tech.JxB_strain, coil_tech.temperature) # A/m^2
            params_Nb3Sn = LTS_scaling(29330000, 28.45, 0.0739, 17.5, -0.7388, -0.5060, -0.0831, 0.8855, 2.169, 2.5,0.0, 1.5, 2.2)
            Jcrit_SC, Bext_Bcrit_ratio = LTS_Jcrit(params_Nb3Sn, Bext, coil_tech.thermal_strain + coil_tech.JxB_strain, coil_tech.temperature) # A/m^2
        elseif coil_tech.material == "NbTi"
            params_NbTi = LTS_scaling(255.3e6, 14.67, -0.002e-2, 8.89, -0.0025, -0.0003, -0.0001, 1.341, 1.555, 2.274, 0.0, 1.758, 2.2) # Table 1, Journal of Phys: Conf. Series, 1559 (2020) 012063
            Jcrit_SC, Bext_Bcrit_ratio = LTS_Jcrit(params_NbTi, Bext, coil_tech.thermal_strain + coil_tech.JxB_strain, coil_tech.temperature) # A/m^2
        elseif coil_tech.material == "KDEMO Nb3Sn"
            Jcrit_SC, Bext_Bcrit_ratio = KDEMO_Nb3Sn_Jcrit(Bext, coil_tech.thermal_strain + coil_tech.JxB_strain, coil_tech.temperature)
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

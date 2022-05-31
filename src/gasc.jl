import JSON
import PeriodicTable: elements

struct GASC
    filename::String
    case::Int
    solution::Dict
    version::Int
end

"""
    GASC(filename::String, case::Int)

Parses GASC output file in json format
"""
function GASC(filename::String, case::Int)
    # from python to julia indexing
    case_j = case + 1
    # parse data
    data = JSON.parsefile(filename)
    # identify version of GASC output
    version = 0
    if "boundaryInnerTF" in keys(data["SOLUTIONS"][case_j]["OUTPUTS"]["numerical profiles"])
        version = 1
    end
    # GASC struct
    gasc = GASC(filename, case, data["SOLUTIONS"][case_j], version)
    # convert list of floats to arrays 
    for item in keys(gasc.solution["OUTPUTS"]["numerical profiles"])
        if item in ["boundaryInnerTF", "boundaryOuterTF"]
        else
            gasc.solution["OUTPUTS"]["numerical profiles"][item] = Vector{Float64}(gasc.solution["OUTPUTS"]["numerical profiles"][item])
        end
    end
    return gasc
end

"""
    case_parameters(gasc::GASC)

Map GASC inputs and solution to FUSE input scalar parameters
"""
function case_parameters(gasc::GASC)
    ini = ParametersInit()
    act = ParametersActor()

    ini.gasc.filename = gasc.filename
    ini.gasc.case = gasc.case
    ini.general.casename = "GASC"
    ini.general.init_from = :scalars

    gasc_2_build(gasc, ini, act)

    gasc_2_equilibrium(gasc, ini, act)

    gasc_2_sources(gasc, ini, act)

    gasc_2_core_profiles(gasc, ini, act)

    return set_new_base!(ini), set_new_base!(act)
end

"""
    gasc_2_core_profiles(gasc::GASC, ini::ParametersInit, act::ParametersActor)

Convert core_profiles information in GASC solution to FUSE `ini` and `act` parameters
"""
function gasc_2_core_profiles(gasc::GASC, ini::ParametersInit, act::ParametersActor)
    gascsol = gasc.solution

    ini.core_profiles.ne_ped = gascsol["OUTPUTS"]["plasma parameters"]["neped"] * 1e20
    ini.core_profiles.greenwald_fraction = gascsol["OUTPUTS"]["plasma parameters"]["greenwaldFraction"]
    ini.core_profiles.T_shaping = 1.8
    i_ped = argmin(abs.(gascsol["OUTPUTS"]["numerical profiles"]["neProf"] .- gascsol["OUTPUTS"]["plasma parameters"]["neped"] / gascsol["OUTPUTS"]["plasma parameters"]["ne0"]))
    ini.core_profiles.w_ped = 1 - gascsol["OUTPUTS"]["numerical profiles"]["rProf"][i_ped]
    ini.core_profiles.zeff = gascsol["OUTPUTS"]["impurities"]["effectiveZ"]
    ini.core_profiles.rot_core = 0.0  # Not in GASC
    ini.core_profiles.bulk = :DT
    ini.core_profiles.helium_fraction = gascsol["INPUTS"]["impurities"]["heliumFraction"]
    ini.core_profiles.impurity = Symbol(elements[Int(gascsol["INPUTS"]["impurities"]["impurityZ"])].symbol)
    ini.core_profiles.ejima = gascsol["INPUTS"]["plasma parameters"]["ejimaCoeff"]
    return ini
end

"""
    gasc_2_equilibrium(gasc::GASC, ini::ParametersInit, act::ParametersActor)

Convert equilibrium information in GASC solution to FUSE `ini` and `act` parameters
"""
function gasc_2_equilibrium(gasc::GASC, ini::ParametersInit, act::ParametersActor)
    gascsol = gasc.solution
    ini.equilibrium.B0 = gascsol["INPUTS"]["conductors"]["magneticFieldOnAxis"]
    ini.equilibrium.R0 = gascsol["INPUTS"]["radial build"]["majorRadius"]
    ini.equilibrium.Z0 = 0.0
    ini.equilibrium.ϵ = 1 / gascsol["INPUTS"]["radial build"]["aspectRatio"]
    ini.equilibrium.κ = gascsol["OUTPUTS"]["plasma parameters"]["elongation"]
    ini.equilibrium.δ = gascsol["INPUTS"]["plasma parameters"]["triangularity"]
    ini.equilibrium.βn = gascsol["OUTPUTS"]["plasma parameters"]["betaN"]
    ini.equilibrium.ip = gascsol["INPUTS"]["plasma parameters"]["plasmaCurrent"] * 1E6
    ini.equilibrium.x_point = gascsol["INPUTS"]["divertor metrics"]["numberDivertors"] > 0
    ini.equilibrium.symmetric = (mod(gascsol["INPUTS"]["divertor metrics"]["numberDivertors"], 2) == 0)
    return ini
end

"""
    gasc_2_sources(gasc::GASC, ini::ParametersInit, act::ParametersActor)

Convert sources (NBI, EC, IC, LH) information in GASC solution to FUSE `ini` and `act` parameters
"""
function gasc_2_sources(gasc::GASC, ini::ParametersInit, act::ParametersActor)
    gascsol = gasc.solution

    inputs = gascsol["INPUTS"]["current drive"]
    outputs = gascsol["OUTPUTS"]["current drive"]

    cd_powers = Float64[]
    ini.nbi.power_launched = Float64[]
    ini.nbi.beam_energy = Float64[]
    ini.nbi.efficiency_conversion = inputs["efficiencyConversionNNBCD"]
    ini.nbi.efficiency_transmission = inputs["efficiencyTransmissionNNBCD"]
    pow = outputs["CDpowerNBCD"] * 1E6 * inputs["NBCDFraction"]
    if pow > 0
        push!(cd_powers, pow)
        push!(ini.nbi.power_launched, pow)
        push!(ini.nbi.beam_energy, 200e3)
    end
    pow = outputs["CDpowerNNBCD"] * 1E6 * inputs["NNBCDFraction"]
    if pow > 0
        push!(cd_powers, pow)
        push!(ini.nbi.power_launched, pow)
        push!(ini.nbi.beam_energy, 1000e3)
    end

    ini.lh_antennas.power_launched = Float64[]
    ini.lh_antennas.efficiency_conversion = inputs["efficiencyConversionLHCD"]
    ini.lh_antennas.efficiency_transmission = inputs["efficiencyTransmissionLHCD"]
    ini.lh_antennas.efficiency_coupling = 1.0 # Not in GASC
    pow = outputs["CDpowerLHCD"] * 1E6 * inputs["LHCDFraction"]
    if pow > 0
        push!(cd_powers, pow)
        push!(ini.lh_antennas.power_launched, pow)
    end
    pow = outputs["CDpowerHICD"] * 1E6 * inputs["HICDFraction"]
    if pow > 0
        push!(cd_powers, pow)
        push!(ini.lh_antennas.power_launched, pow)
    end

    ini.ec_launchers.power_launched = Float64[]
    ini.ec_launchers.efficiency_conversion = inputs["efficiencyConversionECCD"]
    ini.ec_launchers.efficiency_transmission = inputs["efficiencyTransmissionECCD"]
    pow = outputs["CDpowerECCD"] * 1E6 * inputs["ECCDFraction"]
    if pow > 0
        push!(cd_powers, pow)
        push!(ini.ec_launchers.power_launched, pow)
    end

    ini.ic_antennas.power_launched = Float64[]
    ini.ic_antennas.efficiency_conversion = inputs["efficiencyConversionFWCD"]
    ini.ic_antennas.efficiency_transmission = inputs["efficiencyTransmissionFWCD"]
    ini.ic_antennas.efficiency_coupling = 1.0 # Not in GASC
    pow = outputs["CDpowerFWCD"] * 1E6 * inputs["FWCDFraction"]
    if pow > 0
        push!(cd_powers, pow)
        push!(ini.ic_antennas.power_launched, pow)
    end

    # GASC heating power is assumed to be deposited in the core.
    # We use an NBI source as a proxy, which mostly deposits
    # in the core and heats mostly the ions.
    @assert inputs["auxCDPowerFactor"] >= 1.0
    cd_power = sum(cd_powers)
    injected_power = cd_power * inputs["auxCDPowerFactor"]
    heating_power = injected_power - cd_power
    if heating_power > 0
        push!(ini.nbi.power_launched, heating_power)
        push!(ini.nbi.beam_energy, 200e3)
    end
    return ini
end

"""
    gasc_2_build(gasc::GASC, ini::ParametersInit, act::ParametersActor)

Convert radial build information in GASC solution to FUSE `ini` and `act` parameters
"""
function gasc_2_build(gasc::GASC, ini::ParametersInit, act::ParametersActor)
    gascsol = gasc.solution
    ini.build.layers = gasc_2_layers(gascsol)
    ini.build.symmetric = (mod(gascsol["INPUTS"]["divertor metrics"]["numberDivertors"], 2) == 0)

    ini.tf.technology = gasc_2_coil_technology(gasc, :TF)
    ini.oh.technology = gasc_2_coil_technology(gasc, :OH)
    ini.pf_active.technology = gasc_2_coil_technology(gasc, :PF)

    ini.center_stack.bucked = gascsol["INPUTS"]["radial build"]["isBucked"]
    ini.center_stack.noslip = gascsol["INPUTS"]["radial build"]["nonSlip"]
    ini.center_stack.plug = gascsol["INPUTS"]["radial build"]["hasPlug"]

    ini.oh.flattop_duration = gascsol["INPUTS"]["plasma parameters"]["flattopDuration"]
    return ini
end

"""
    function gasc_2_layers(gascsol::Dict)

Convert GASC ["OUTPUTS"]["radial build"] to FUSE build layers dictionary
"""
function gasc_2_layers(gascsol::Dict)
    gascrb = gascsol["OUTPUTS"]["radial build"]

    layers = DataStructures.OrderedDict()
    mapper = Dict(
        "OH" => "OH",
        "TF" => "TF",
        "LTShield" => "low_temp_shield",
        "VV" => "vacuum_vessel",
        "HTShield" => "high_temp_shield",
        "Blanket" => "blanket",
        "Plasma" => "plasma")

    gasc_layers = [
        "RiOH",
        "RoOH",
        "RiInnerTF",
        "RoInnerTF",
        "RiInnerLTShield",
        "RoInnerLTShield",
        "RiInnerVV",
        "RoInnerVV",
        "RiInnerHTShield",
        "RoInnerHTShield",
        "RiInnerBlanket",
        "RoInnerBlanket",
        "RiPlasma",
        "RoPlasma",
        "RiOuterBlanket",
        "RoOuterBlanket",
        "RiOuterHTShield",
        "RoOuterHTShield",
        "RiOuterVV",
        "RoOuterVV",
        "RiOuterLTShield",
        "RoOuterLTShield",
        "RiOuterTF",
        "RoOuterTF"]

    layers["gap_OH"] = gascrb["RiOH"]
    for k in 2:length(gasc_layers)
        g1 = gasc_layers[k-1]
        g2 = gasc_layers[k]
        d = gascrb[replace(g2, "InnerTF" => "TF")] - gascrb[replace(g1, "InnerTF" => "TF")]
        f1 = mapper[replace(replace(replace(replace(g1, "Ri" => ""), "Ro" => ""), "Inner" => ""), "Outer" => "")]
        f2 = mapper[replace(replace(replace(replace(g2, "Ri" => ""), "Ro" => ""), "Inner" => ""), "Outer" => "")]
        if startswith(g1, "Ri")
            if contains(g1, "Inner")
                f1 = "hfs_" * f1
            elseif contains(g1, "Outer")
                f1 = "lfs_" * f1
            end
            f = f1
            f = replace(f, r".fs_gap_TF_OH" => "gap_TF_OH")
            layers[f] = d
        elseif startswith(g1, "Ro") && d > 0
            if contains(g2, "Outer")
                f = "lfs_gap_$(f1)_$(f2)"
            else
                f = "hfs_gap_$(f2)_$(f1)"
            end
            f = replace(f, r".fs_gap_TF_OH" => "gap_TF_OH")
            layers[f] = d
        end
    end
    layers["gap_cryostat"] = layers["OH"] * 4
    layers["cryostat"] = layers["lfs_low_temp_shield"] / 2.0

    for k in collect(keys(layers))
        if contains(k, "_gap_") && contains(k, "_plasma_")
            layers["plasma"] += layers[k]
            delete!(layers, k)
        end
    end

    if false # to remove gap that prevents bucking
        for k in collect(keys(layers))
            if k == "gap_TF_OH"
                layers["OH"] += layers["gap_TF_OH"]
                delete!(layers, k)
            end
        end
    end

    return layers
end

"""
    gasc_2_coil_technology(gasc::GASC, coil_type::Symbol)

Convert coil technology information in GASC solution to FUSE `coil_technology` data structure
"""
function gasc_2_coil_technology(gasc::GASC, coil_type::Symbol)
    gascsol = gasc.solution
    if coil_type ∉ [:OH, :TF, :PF]
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
        if coil_type == :PF # assume PF coils are always LTS
            coil_tech = coil_technology(:LTS)
        end
        if "thermalStrain$coil_type" ∉ keys(gascsol["INPUTS"]["conductors"])
            coil_tech.thermal_strain = 0.0
            coil_tech.JxB_strain = 0.0
        else
            coil_tech.thermal_strain = gascsol["INPUTS"]["conductors"]["thermalStrain$coil_type"]
            coil_tech.JxB_strain = gascsol["INPUTS"]["conductors"]["structuralStrain$coil_type"]
        end
    end
    if "fractionVoid$coil_type" ∉ keys(gascsol["INPUTS"]["conductors"])
        coil_tech.fraction_void = gascsol["INPUTS"]["conductors"]["fractionVoidOH"]
        coil_tech.fraction_stainless = gascsol["INPUTS"]["conductors"]["fractionStainlessOH"]
    else
        coil_tech.fraction_void = gascsol["INPUTS"]["conductors"]["fractionVoid$coil_type"]
        coil_tech.fraction_stainless = gascsol["INPUTS"]["conductors"]["fractionStainless$coil_type"]
    end
    return set_new_base!(coil_tech)
end

function compare(dd::IMAS.dd, gasc::GASC)
    df = DataFrames.DataFrame(code=["FUSE", "GASC"])

    # collisionless bootstrap coefficient
    FUSE = IMAS.collisionless_bootstrap_coefficient(dd)
    GASC = gasc.solution["INPUTS"]["plasma parameters"]["user_bootstrapCoefficient"]
    df.Cbs = [FUSE, GASC]

    # fusion power [MW]
    FUSE = IMAS.fusion_power(dd.core_profiles.profiles_1d[]) / 1E6
    GASC = gasc.solution["OUTPUTS"]["power balance"]["powerFusion"]
    df.Pfusion = [FUSE, GASC]

    df
end
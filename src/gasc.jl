import JSON

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
    ini = InitParameters()
    act = ActorParameters()

    ini.gasc.filename = gasc.filename
    ini.gasc.case = gasc.case
    ini.general.casename = "GASC"
    ini.general.init_from = :gasc

    gasc_2_build(gasc, ini, act)

    gasc_2_equilibrium(gasc, ini, act)

    gasc_2_sources(gasc, ini, act)

    return set_new_base!(ini), set_new_base!(act)
end

function gasc_2_equilibrium(gasc::GASC, ini::InitParameters, act::ActorParameters)
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

function gasc_2_sources(gasc::GASC, ini::InitParameters, act::ActorParameters)
    gascsol = gasc.solution

    inputs = gascsol["INPUTS"]["current drive"]
    outputs = gascsol["OUTPUTS"]["current drive"]

    cd_powers = Float64[]
    ini.nbi.power_launched = Float64[]
    ini.nbi.beam_energy = Float64[]
    if outputs["CDpowerNBCD"] > 0
        pow = outputs["CDpowerNBCD"] * 1E6 * inputs["NBCDFraction"]
        push!(cd_powers, pow)
        push!(ini.nbi.power_launched, pow)
        push!(ini.nbi.beam_energy, 200e3)
    end
    if outputs["CDpowerNNBCD"] > 0
        pow = outputs["CDpowerNNBCD"] * 1E6 * inputs["NNBCDFraction"]
        push!(cd_powers, pow)
        push!(ini.nbi.power_launched, pow)
        push!(ini.nbi.beam_energy, 1000e3)
    end

    ini.lh.power_launched = Float64[]
    if outputs["CDpowerLHCD"] > 0
        pow = outputs["CDpowerLHCD"] * 1E6 * inputs["LHCDFraction"]
        push!(cd_powers, pow)
        push!(ini.lh.power_launched, pow)
    end
    if outputs["CDpowerHICD"] > 0
        pow = outputs["CDpowerHICD"] * 1E6 * inputs["HICDFraction"]
        push!(cd_powers, pow)
        push!(ini.lh.power_launched, pow)
    end

    ini.ec.power_launched = Float64[]
    if outputs["CDpowerECCD"] > 0
        pow = outputs["CDpowerECCD"] * 1E6 * inputs["ECCDFraction"]
        push!(cd_powers, pow)
        push!(ini.ec.power_launched, pow)
    end

    ini.ic.power_launched = Float64[]
    if outputs["CDpowerFWCD"] > 0
        pow = outputs["CDpowerFWCD"] * 1E6 * inputs["FWCDFraction"]
        push!(cd_powers, pow)
        push!(ini.ic.power_launched, pow)
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

function gasc_2_build(gasc::GASC, ini::InitParameters, act::ActorParameters)
    gascsol = gasc.solution
    ini.build.layers = gasc_to_layers(gascsol)
    ini.build.symmetric = (mod(gascsol["INPUTS"]["divertor metrics"]["numberDivertors"], 2) == 0)

    ini.tf.technology = coil_technology(gasc, :TF)
    ini.oh.technology = coil_technology(gasc, :OH)
    ini.pf_active.technology = coil_technology(gasc, :PF)

    ini.center_stack.bucked = gascsol["INPUTS"]["radial build"]["isBucked"]
    ini.center_stack.noslip = gascsol["INPUTS"]["radial build"]["nonSlip"]
    ini.center_stack.plug = gascsol["INPUTS"]["radial build"]["hasPlug"]

    ini.oh.flattop_duration = gascsol["INPUTS"]["plasma parameters"]["flattopDuration"]
    return ini
end

"""
    function gasc_to_layers(gascsol::Dict)

Convert GASC ["INPUTS"]["radial build"] to FUSE build layers dictionary
"""
function gasc_to_layers(gascsol::Dict)
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
            layers[f1] = d
        elseif startswith(g1, "Ro") && d > 0
            if contains(g2, "Outer")
                f = "lfs_gap_$(f1)_$(f2)"
            else
                f = "hfs_gap_$(f2)_$(f1)"
            end
            layers[f] = d
        end
    end
    if "hfs_gap_TF_OH" in keys(layers)
        layers["lfs_gap_TF_OH"] = layers["hfs_gap_TF_OH"]
    end
    layers["gap_cryostat"] = layers["gap_OH"] * 3

    return layers
end

"""
    coil_technology(gasc::GASC, coil_type::Symbol)

Return coil parameters from GASC solution and coil type [:OH, :TF, :PF]
"""
function coil_technology(gasc::GASC, coil_type::Symbol)
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
        if "thermalStrain$coil_type" ∉ keys(gascsol["INPUTS"]["conductors"])
            coil_tech.thermal_strain = 0.0
            coil_tech.JxB_strain = 0.0
        else
            coil_tech.thermal_strain = gascsol["INPUTS"]["conductors"]["thermalStrain$coil_type"]
            coil_tech.JxB_strain = gascsol["INPUTS"]["conductors"]["structuralStrain$coil_type"]
        end
    end
    if "fractionVoid$coil_type" ∉ keys(gascsol["INPUTS"]["conductors"])
        coil_tech.fraction_void = 0.1
        coil_tech.fraction_stainless = 0.5
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
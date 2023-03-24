import JSON
import PeriodicTable: elements

mutable struct GASC
    filename::String
    data::Dict
    case::Int
    version::Int
end

function Base.getproperty(gasc::GASC, field::Symbol)
    if field == :solution
        return getfield(gasc, :data)["SOLUTIONS"][gasc.case+1]
    elseif field == :outputs
        return getfield(gasc, :data)["SOLUTIONS"][gasc.case+1]["OUTPUTS"]
    elseif field == :inputs
        return getfield(gasc, :data)["SOLUTIONS"][gasc.case+1]["INPUTS"]
    elseif field == :constraints
        return getfield(gasc, :data)["SETUP"]["SETTINGS"]["constraints"]
    else
        return getfield(gasc, field)
    end
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
    gasc = GASC(filename, data, case, version)
    # convert list of floats to arrays
    for item in keys(gasc.outputs["numerical profiles"])
        if item ∉ ["boundaryInnerTF", "boundaryOuterTF"]
            gasc.outputs["numerical profiles"][item] = Vector{Float64}(gasc.outputs["numerical profiles"][item])
        end
    end
    return gasc
end

"""
    case_parameters(gasc::GASC)

Map GASC inputs and solution to FUSE input scalar parameters
"""
function case_parameters(gasc::GASC)
    ini = ParametersInits()
    act = ParametersActors()

    ini.gasc.filename = gasc.filename
    ini.gasc.case = gasc.case
    ini.general.casename = gasc.data["LABEL"]
    ini.general.init_from = :scalars

    gasc_2_build(gasc, ini, act)

    gasc_2_target(gasc, ini, act)

    gasc_2_equilibrium(gasc, ini, act)

    gasc_2_sources(gasc, ini, act)

    gasc_2_core_profiles(gasc, ini, act)

    return set_new_base!(ini), set_new_base!(act)
end

"""
    gasc_2_core_profiles(gasc::GASC, ini::ParametersAllInits, act::ParametersAllActors)

Convert core_profiles information in GASC solution to FUSE `ini` and `act` parameters
"""
function gasc_2_core_profiles(gasc::GASC, ini::ParametersAllInits, act::ParametersAllActors)
    ini.core_profiles.greenwald_fraction = gasc.outputs["plasma parameters"]["greenwaldFraction"]
    # ini.core_profiles.ne_ped = gasc.outputs["plasma parameters"]["neped"] * 1e20
    # i_ped = argmin(abs.(gasc.outputs["numerical profiles"]["neProf"] .- gasc.outputs["plasma parameters"]["neped"] / gasc.outputs["plasma parameters"]["ne0"]))
    # ini.core_profiles.w_ped = 1 - gasc.outputs["numerical profiles"]["rProf"][i_ped]
    ini.core_profiles.T_shaping = 1.8
    ini.core_profiles.zeff = gasc.outputs["impurities"]["effectiveZ"]
    ini.core_profiles.rot_core = 0.0  # Not in GASC
    ini.core_profiles.bulk = :DT
    ini.core_profiles.helium_fraction = gasc.inputs["impurities"]["heliumFraction"]
    ini.core_profiles.impurity = Symbol(elements[Int(gasc.inputs["impurities"]["impurityZ"])].symbol)
    ini.core_profiles.ejima = gasc.inputs["plasma parameters"]["ejimaCoeff"]
    return ini
end

"""
    gasc_2_equilibrium(gasc::GASC, ini::ParametersAllInits, act::ParametersAllActors)

Convert equilibrium information in GASC solution to FUSE `ini` and `act` parameters
"""
function gasc_2_equilibrium(gasc::GASC, ini::ParametersAllInits, act::ParametersAllActors)
    ini.equilibrium.boundary_from = :scalars
    ini.equilibrium.B0 = gasc.inputs["conductors"]["magneticFieldOnAxis"]
    ini.equilibrium.R0 = gasc.inputs["radial build"]["majorRadius"]
    ini.equilibrium.Z0 = 0.0
    ini.equilibrium.ϵ = 1.0 / gasc.inputs["radial build"]["aspectRatio"]
    ini.equilibrium.κ = gasc.outputs["plasma parameters"]["elongation"]
    ini.equilibrium.δ = gasc.inputs["plasma parameters"]["triangularity"]

    Pavg = gasc.outputs["plasma parameters"]["pressureVolAvg"]
    V = gasc.outputs["plasma parameters"]["plasmaVolume"]
    vol = gasc.outputs["numerical profiles"]["volumeProf"] .* V
    P1 = sum(IMAS.gradient(vol) .* LinRange(1.0, 0.0, length(vol))) / V
    ini.equilibrium.pressure_core = Pavg / P1

    ini.equilibrium.ip = gasc.inputs["plasma parameters"]["plasmaCurrent"] * 1E6
    ini.equilibrium.xpoints_number = Int(round(gasc.inputs["divertor metrics"]["numberDivertors"]))
    return ini
end

"""
    gasc_2_sources(gasc::GASC, ini::ParametersAllInits, act::ParametersAllActors)

Convert sources (NBI, EC, IC, LH) information in GASC solution to FUSE `ini` and `act` parameters
"""
function gasc_2_sources(gasc::GASC, ini::ParametersAllInits, act::ParametersAllActors)
    inputs = gasc.inputs["current drive"]
    outputs = gasc.outputs["current drive"]

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
    gasc_2_build(gasc::GASC, ini::ParametersAllInits, act::ParametersAllActors)

Convert radial build information in GASC solution to FUSE `ini` and `act` parameters
"""
function gasc_2_build(gasc::GASC, ini::ParametersAllInits, act::ParametersAllActors)
    ini.build.layers = gasc_2_layers(gasc)
    ini.build.symmetric = (mod(gasc.inputs["divertor metrics"]["numberDivertors"], 2) == 0)

    ini.tf.technology = gasc_2_coil_technology(gasc, :TF)
    ini.oh.technology = gasc_2_coil_technology(gasc, :OH)
    ini.pf_active.technology = gasc_2_coil_technology(gasc, :PF)

    ini.center_stack.bucked = gasc.inputs["radial build"]["isBucked"]
    ini.center_stack.noslip = gasc.inputs["radial build"]["nonSlip"]
    ini.center_stack.plug = gasc.inputs["radial build"]["hasPlug"]
    return ini
end

"""
    gasc_2_target(gasc::GASC, ini::ParametersAllInits, act::ParametersAllActors)

Convert nominal target design information in GASC solution to FUSE `ini` and `act` parameters
"""
function gasc_2_target(gasc::GASC, ini::ParametersAllInits, act::ParametersAllActors)
    ini.requirements.flattop_duration = float(gasc.inputs["plasma parameters"]["flattopDuration"])
    ini.requirements.power_electric_net = gasc.constraints["lowerOutputConstraints"]["powerNet"] * 1E6
    return ini
end

"""
    function gasc_2_layers(gasc::GASC)

Convert GASC ["OUTPUTS"]["radial build"] to FUSE build layers dictionary
"""
function gasc_2_layers(gasc::GASC)
    gascrb = gasc.outputs["radial build"]

    layers = OrderedCollections.OrderedDict()
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
        if startswith(g1, "Ri") # hfs
            if contains(g1, "Inner")
                f1 = "hfs_" * f1
            elseif contains(g1, "Outer")
                f1 = "lfs_" * f1
            end
            f = f1
            f = replace(f, r".fs_gap_TF_OH" => "gap_TF_OH")
            layers[f] = d

        elseif startswith(g1, "Ro") && (d > 0) # lfs
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

    return layers
end

"""
    gasc_buck_OH_TF!(layers::OrderedCollections.OrderedDict)

Remove gap between OH and TF to allow bucking (gap gets added to OH thickness)
"""
function gasc_buck_OH_TF!(layers::OrderedCollections.OrderedDict)
    for k in collect(keys(layers))
        if k == "gap_TF_OH"
            layers["OH"] += layers["gap_TF_OH"]
            delete!(layers, k)
        end
    end
    return layers
end

"""
    gasc_add_wall_layers!(layers::OrderedCollections.OrderedDict, thickness::Float64)

Add wall layer of given thickness expressed [meters] (gets subtracted from blanket layer)
"""
function gasc_add_wall_layers!(layers::OrderedCollections.OrderedDict; thickness::Float64)
    tmp = OrderedCollections.OrderedDict()
    for layer in keys(layers)
        if layer == "hfs_blanket"
            tmp[layer] = layers[layer] - thickness
            tmp["hfs_first_wall"] = thickness
        elseif layer == "lfs_blanket"
            tmp["lfs_first_wall"] = thickness
            tmp[layer] = layers[layer] - thickness
        elseif layer == "hfs_vacuum_vessel"
            tmp["hfs_vacuum_vessel_wall_outer"] = thickness
            tmp[layer] = layers[layer] - 2 * thickness
            tmp["hfs_vacuum_vessel_wall_inner"] = thickness
        elseif layer == "lfs_vacuum_vessel"
            tmp["lfs_vacuum_vessel_wall_inner"] = thickness
            tmp[layer] = layers[layer] - 2 * thickness
            tmp["lfs_vacuum_vessel_wall_outer"] = thickness
        else
            tmp[layer] = layers[layer]
        end
    end
    empty!(layers)
    for layer in keys(tmp)
        layers[layer] = tmp[layer]
    end
    return layers
end

"""
    gasc_2_coil_technology(gasc::GASC, coil_type::Symbol)

Convert coil technology information in GASC solution to FUSE `coil_technology` data structure
"""
function gasc_2_coil_technology(gasc::GASC, coil_type::Symbol)
    if coil_type ∉ [:OH, :TF, :PF]
        error("Supported coil type are [:OH, :TF, :PF]")
    end
    if gasc.inputs["conductors"]["superConducting"] == "copper"
        coil_tech = coil_technology(:copper)
    else
        if gasc.inputs["conductors"]["superConducting"] == "LTS"
            coil_tech = coil_technology(:LTS)
        elseif gasc.inputs["conductors"]["superConducting"] == "HTS"
            coil_tech = coil_technology(:HTS)
        end
        if coil_type == :PF # assume PF coils are always LTS
            coil_tech = coil_technology(:LTS)
        end
        if "thermalStrain$coil_type" ∉ keys(gasc.inputs["conductors"])
            coil_tech.thermal_strain = 0.0
            coil_tech.JxB_strain = 0.0
        else
            coil_tech.thermal_strain = float(gasc.inputs["conductors"]["thermalStrain$coil_type"])
            coil_tech.JxB_strain = float(gasc.inputs["conductors"]["structuralStrain$coil_type"])
        end
    end
    if "fractionVoid$coil_type" ∉ keys(gasc.inputs["conductors"])
        coil_tech.fraction_void = float(gasc.inputs["conductors"]["fractionVoidOH"])
        coil_tech.fraction_stainless = float(gasc.inputs["conductors"]["fractionStainlessOH"])
    else
        coil_tech.fraction_void = float(gasc.inputs["conductors"]["fractionVoid$coil_type"])
        coil_tech.fraction_stainless = float(gasc.inputs["conductors"]["fractionStainless$coil_type"])
    end
    return set_new_base!(coil_tech)
end

function compare(dd::IMAS.dd, gasc::GASC)
    df = DataFrames.DataFrame(code=["FUSE", "GASC"])

    # collisionless bootstrap coefficient
    FUSE = IMAS.collisionless_bootstrap_coefficient(dd)
    GASC = gasc.inputs["plasma parameters"]["user_bootstrapCoefficient"]
    df.Cbs = [FUSE, GASC]

    # fusion power [MW]
    FUSE = IMAS.fusion_power(dd.core_profiles) / 1E6
    GASC = gasc.outputs["power balance"]["powerFusion"]
    df.Pfusion = [FUSE, GASC]

    df
end
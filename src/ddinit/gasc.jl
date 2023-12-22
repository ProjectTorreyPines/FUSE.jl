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
    data = JSON.parsefile(replace(filename, r"__FUSE__" => __FUSE__))
    # identify version of GASC output
    version = 0
    if "boundaryInnerTF" in keys(data["SOLUTIONS"][case_j]["OUTPUTS"]["numerical profiles"])
        version = 1
    end
    # GASC struct
    gasc = GASC(filename, data, case, version)
    # convert list of floats to arrays
    for item in keys(gasc.outputs["numerical profiles"])
        if item ∉ ("boundaryInnerTF", "boundaryOuterTF")
            gasc.outputs["numerical profiles"][item] = Vector{Float64}(gasc.outputs["numerical profiles"][item])
        end
    end
    return gasc
end

"""
    case_parameters(gasc::GASC)

Map GASC inputs and solution to FUSE input scalar parameters
"""
function case_parameters(gasc::GASC; add_wall_layers::Float64=0.0)
    ini = ParametersInits()
    act = ParametersActors()

    ini.gasc.filename = gasc.filename
    ini.gasc.case = gasc.case
    ini.general.casename = gasc.data["LABEL"]
    ini.general.init_from = :scalars

    gasc_2_build(gasc, ini, act; add_wall_layers)

    gasc_2_requirement(gasc, ini, act)

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
    ini.core_profiles.T_ratio = 1.0
    ini.core_profiles.T_shaping = 1.8
    ini.core_profiles.n_shaping = 0.9
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
    P1 = sum(IMAS.gradient(vol) .* range(1.0, 0.0, length(vol))) / V
    ini.equilibrium.pressure_core = Pavg / P1

    ini.equilibrium.ip = gasc.inputs["plasma parameters"]["plasmaCurrent"] * 1E6
    nx = Int(gasc.inputs["divertor metrics"]["numberDivertors"])
    if nx == 0
        ini.equilibrium.xpoints = :none
    elseif nx == 1
        ini.equilibrium.xpoints = :lower
    elseif nx == 2
        ini.equilibrium.xpoints = :double
    else
        error("Invalid GASC numberDivertors")
    end

    return ini
end

"""
    gasc_2_sources(gasc::GASC, ini::ParametersAllInits, act::ParametersAllActors)

Convert sources (NB, EC, IC, LH) information in GASC solution to FUSE `ini` and `act` parameters
"""
function gasc_2_sources(gasc::GASC, ini::ParametersAllInits, act::ParametersAllActors)
    inputs = gasc.inputs["current drive"]
    outputs = gasc.outputs["current drive"]

    T = Float64
    cd_powers = T[]

    push!(ini.nb_unit, )

    for (energy,pow) in ((200e3, outputs["CDpowerNBCD"] * 1E6 * inputs["NBCDFraction"]), (1000e3, outputs["CDpowerNNBCD"] * 1E6 * inputs["NNBCDFraction"]))
        if pow > 0.0
            nb_unit = FUSEparameters__nb_unit{T}()
            nb_unit.efficiency_conversion = inputs["efficiencyConversionNNBCD"]
            nb_unit.efficiency_transmission = inputs["efficiencyTransmissionNNBCD"]
            nb_unit.beam_energy = energy
            nb_unit.power_launched = pow

            push!(ini.nb_unit, nb_unit)
            push!(cd_powers, pow)
        end
    end

    pow = outputs["CDpowerECCD"] * 1E6 * inputs["ECCDFraction"]
    if pow > 0.0
        ec_launcher = FUSEparameters__ec_launcher{T}()
        ec_launcher.efficiency_conversion = inputs["efficiencyConversionECCD"]
        ec_launcher.efficiency_transmission = inputs["efficiencyTransmissionECCD"]
        ec_launcher.power_launched = pow

        push!(ini.ec_launcher, ec_launcher)
        push!(cd_powers, pow)
    end

    pow = outputs["CDpowerFWCD"] * 1E6 * inputs["FWCDFraction"]
    if pow > 0.0
        ic_antenna = FUSEparameters__ic_antenna{T}()
        ic_antenna.efficiency_conversion = inputs["efficiencyConversionFWCD"]
        ic_antenna.efficiency_transmission = inputs["efficiencyTransmissionFWCD"]
        ic_antenna.efficiency_coupling = 1.0 # Not in GASC
        ic_antenna.power_launched = pow

        push!(ini.ic_antenna, ic_antenna)
        push!(cd_powers, pow)
    end


    for pow in (outputs["CDpowerLHCD"] * 1E6 * inputs["LHCDFraction"], outputs["CDpowerHICD"] * 1E6 * inputs["HICDFraction"])
        if pow > 0.0
            lh_antenna = FUSEparameters__lh_antenna{T}()
            lh_antenna.efficiency_conversion = inputs["efficiencyConversionLHCD"]
            lh_antenna.efficiency_transmission = inputs["efficiencyTransmissionLHCD"]
            lh_antenna.efficiency_coupling = 1.0 # Not in GASC
            lh_antenna.power_launched = pow

            push!(ini.lh_antenna, lh_antenna)
            push!(cd_powers, pow)
        end
    end

    # GASC heating power is assumed to be deposited in the core.
    # We use an NBI source as a proxy, which mostly deposits
    # in the core and heats mostly the ions.
    @assert inputs["auxCDPowerFactor"] >= 1.0
    cd_power = sum(cd_powers)
    injected_power = cd_power * inputs["auxCDPowerFactor"]
    heating_power = injected_power - cd_power
    if heating_power > 0.0
        nb_unit = FUSEparameters__nb_unit{T}()
        nb_unit.efficiency_conversion = inputs["efficiencyConversionNNBCD"]
        nb_unit.efficiency_transmission = inputs["efficiencyTransmissionNNBCD"]
        nb_unit.beam_energy = 200e3
        nb_unit.power_launched = heating_power

        push!(ini.nb_unit, nb_unit)
    end

    # # set power_launched to missing or scalar if number of actuators is 0 or 1
    # for heating_scheme in (ini.nb_unit, ini.ec_launcher, ini.ic_antenna, ini.lh_antenna)
    #     if isempty(heating_scheme.power_launched)
    #         heating_scheme.power_launched = missing
    #     elseif length(heating_scheme.power_launched) == 1
    #         heating_scheme.power_launched = heating_scheme.power_launched[1]
    #     end
    # end

    return ini
end

"""
    gasc_2_build(gasc::GASC, ini::ParametersAllInits, act::ParametersAllActors)

Convert radial build information in GASC solution to FUSE `ini` and `act` parameters
"""
function gasc_2_build(gasc::GASC, ini::ParametersAllInits, act::ParametersAllActors; add_wall_layers::Float64)
    layers = gasc_2_layers(gasc)
    if add_wall_layers > 0.0
        gasc_add_wall_layers!(layers; thickness=add_wall_layers)
        ini.build.n_first_wall_conformal_layers = 2
    end
    ini.build.layers = layers

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
    gasc_2_requirement(gasc::GASC, ini::ParametersAllInits, act::ParametersAllActors)

Convert requirements of GASC solution to FUSE `ini` and `act` parameters
"""
function gasc_2_requirement(gasc::GASC, ini::ParametersAllInits, act::ParametersAllActors)
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

    layers = OrderedCollections.OrderedDict{String,Float64}()
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

    # convert keys from string to symbols
    sym_layers = OrderedCollections.OrderedDict{Symbol,Float64}()
    for k in keys(layers)
        sym_layers[Symbol(k)] = layers[k]
    end

    return sym_layers
end

"""
    gasc_buck_OH_TF!(layers::OrderedCollections.OrderedDict{Symbol,Float64})

Remove gap between OH and TF to allow bucking (gap gets added to OH thickness)
"""
function gasc_buck_OH_TF!(layers::OrderedCollections.OrderedDict{Symbol,Float64})
    for layer in collect(keys(layers))
        if layer == :gap_TF_OH
            layers[:OH] += layers[:gap_TF_OH]
            delete!(layers, k)
        end
    end
    return layers
end

"""
    gasc_add_wall_layers!(layers::OrderedCollections.OrderedDict{Symbol,Float64}; thickness::Float64)

Add wall layer of given thickness expressed [meters] (gets subtracted from blanket layer)
"""
function gasc_add_wall_layers!(layers::OrderedCollections.OrderedDict{Symbol,Float64}; thickness::Float64)
    tmp = OrderedCollections.OrderedDict{Symbol,Float64}()
    for layer in keys(layers)
        if layer == :hfs_blanket
            tmp[layer] = layers[layer] - thickness
            tmp[:hfs_first_wall] = thickness
        elseif layer == :lfs_blanket
            tmp[:lfs_first_wall] = thickness
            tmp[layer] = layers[layer] - thickness
        elseif layer == :hfs_vacuum_vessel
            tmp[:hfs_vacuum_vessel_wall_outer] = thickness
            tmp[layer] = layers[layer] - 2 * thickness
            tmp[:hfs_vacuum_vessel_wall_inner] = thickness
        elseif layer == :lfs_vacuum_vessel
            tmp[:lfs_vacuum_vessel_wall_inner] = thickness
            tmp[layer] = layers[layer] - 2 * thickness
            tmp[:lfs_vacuum_vessel_wall_outer] = thickness
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

Convert coil technology information in GASC solution to FUSE `coil_technology` Symbol
"""
function gasc_2_coil_technology(gasc::GASC, coil_type::Symbol)
    if coil_type ∉ (:OH, :TF, :PF)
        error("Supported coil type are [:OH, :TF, :PF]")
    end
    if gasc.inputs["conductors"]["superConducting"] == "copper"
        coil_tech = :copper
    else
        if gasc.inputs["conductors"]["superConducting"] == "LTS"
            coil_tech = :Nb3Sn
        elseif gasc.inputs["conductors"]["superConducting"] == "HTS"
            coil_tech = :HTS
        end
        if coil_type == :PF # assume PF coils are always LTS
            coil_tech = :Nb3Sn
        end
    end
    return coil_tech
end

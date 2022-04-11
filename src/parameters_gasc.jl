"""
    case_parameters(gasc::GASC)

Map GASC inputs and solution to FUSE input scalar parameters
"""
function InitParameters(gasc::GASC; no_small_gaps::Bool=true, vacuum_vessel::Float64=0.1)
    ini = InitParameters()

    ini.gasc.filename = gasc.filename
    ini.gasc.case = gasc.case
    ini.general.casename = gasc.solution["INPUTS"]["NAME"]["device_name"]
    ini.general.init_from = :gasc

    gasc_2_build(ini, gasc; no_small_gaps, vacuum_vessel)

    gasc_2_equilibrium(ini, gasc)

    gasc_2_sources_par(ini, gasc)

    return set_new_base!(ini)
end

function gasc_2_equilibrium(ini::InitParameters, gasc::GASC)
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

function gasc_2_sources_par(ini::InitParameters, gasc::GASC)
    gascsol = gasc.solution

    injected_power = gascsol["OUTPUTS"]["current drive"]["powerAux"] * 1E6
    @assert gascsol["INPUTS"]["current drive"]["auxCDPowerFactor"] >= 1.0
    cd_power = injected_power / gascsol["INPUTS"]["current drive"]["auxCDPowerFactor"]
    heating_power = injected_power - cd_power
    plug_power = injected_power / gascsol["INPUTS"]["power efficiency"]["efficiencyAux"]

    cd_powers = Dict()
    for system in ["NNB", "NB", "LH", "FW", "EC", "HI"]
        cd_powers[system] = cd_power * gascsol["INPUTS"]["current drive"]["$(system)CDFraction"]
    end

    ini.nbi.power_launched = Float64[]
    ini.nbi.beam_energy = Float64[]
    if heating_power >0
        push!(ini.nbi.power_launched, heating_power)
        push!(ini.nbi.beam_energy, 200e3)
    end
    if cd_powers["NB"] >0
        push!(ini.nbi.power_launched, cd_powers["NB"])
        push!(ini.nbi.beam_energy, 200e3)
    end
    if cd_powers["NNB"] >0
        push!(ini.nbi.power_launched, cd_powers["NNB"])
        push!(ini.nbi.beam_energy, 1000e3)
    end
    ini.lh.power_launched = Float64[]
    if cd_powers["LH"] > 0
        push!(ini.lh.power_launched, cd_powers["LH"])
    end
    if cd_powers["HI"] > 0
        push!(ini.lh.power_launched, cd_powers["HI"])
    end
    ini.ic.power_launched = cd_powers["FW"]
    ini.ec.power_launched = cd_powers["EC"]

    return ini
end

function gasc_2_build(ini::InitParameters, gasc::GASC; no_small_gaps::Bool, vacuum_vessel::Float64)
    gascsol = gasc.solution
    ini.build.layers = gasc_to_layers(gascsol["INPUTS"]["radial build"]; no_small_gaps, vacuum_vessel)
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
    function gasc_to_layers(
        gascrb::Dict;
        no_small_gaps::Bool = true,
        vacuum_vessel::Float64 = 0.1)

Convert GASC ["INPUTS"]["radial build"] to FUSE build layers dictionary
"""
function gasc_to_layers(
    gascrb::Dict;
    no_small_gaps::Bool,
    vacuum_vessel::Float64)

    majorRadius = gascrb["majorRadius"]
    aspectRatio = gascrb["aspectRatio"]
    minorRadius = majorRadius / aspectRatio
    innerPlasmaRadius = majorRadius - minorRadius
    norm = innerPlasmaRadius

    layers = DataStructures.OrderedDict()
    for run in 1:2
        layers["OH"] = gascrb["rbOH"] * norm

        layers["hfs_gap_TF"] = gascrb["gapTFOH"] * norm
        layers["hfs_TF"] = gascrb["rbTF"] * norm
        if no_small_gaps
            layers["hfs_TF"] += layers["hfs_gap_TF"]
            pop!(layers, "hfs_gap_TF")
        end

        if vacuum_vessel > 0.0
            layers["gap_hfs_vacuum_vessel"] = gascrb["rbInnerBlanket"] * norm * vacuum_vessel
        end

        layers["hfs_gap_shield"] = gascrb["gapBlanketCoil"] * norm
        layers["hfs_shield"] = gascrb["rbInnerShield"] * norm
        if no_small_gaps
            layers["hfs_shield"] += layers["hfs_gap_shield"]
            pop!(layers, "hfs_gap_shield")
        end
        layers["hfs_blanket"] = gascrb["rbInnerBlanket"] * norm * (1 - vacuum_vessel)

        layers["hfs_wall"] = gascrb["gapInnerBlanketWall"] * norm

        if run == 1
            between_gapOH_and_plasma = sum(values(layers))
            empty!(layers)
            layers["gap_OH"] = innerPlasmaRadius - between_gapOH_and_plasma
        end
    end

    layers["plasma"] = (gascrb["majorRadius"] - sum(values(layers))) * 2 + (gascrb["gapIn"] + gascrb["gapOut"]) * norm
    layers["lfs_wall"] = gascrb["gapOuterBlanketWall"] * norm

    layers["lfs_blanket"] = gascrb["rbOuterBlanket"] * norm * (1 - vacuum_vessel)
    layers["lfs_shield"] = gascrb["rbOuterShield"] * norm
    layers["lfs_gap_shield"] = gascrb["gapBlanketCoil"] * norm
    if no_small_gaps
        layers["lfs_shield"] += layers["lfs_gap_shield"]
        pop!(layers, "lfs_gap_shield")
    end

    if vacuum_vessel > 0.0
        layers["gap_lfs_vacuum_vessel"] = gascrb["rbOuterBlanket"] * norm * vacuum_vessel
    end

    layers["lfs_TF"] = layers["hfs_TF"]
    layers["lfs_gap_TF"] = gascrb["gapTFOH"] * norm
    if no_small_gaps
        layers["lfs_TF"] += layers["lfs_gap_TF"]
        pop!(layers, "lfs_gap_TF")
    end

    layers["gap_cryostat"] = layers["gap_OH"] * 3

    # thin layers can cause LibGEOS to crash
    # if wall is too thin, then thicken it at the expense of the blanket
    min_fraction_thin_wall = 0.02
    if no_small_gaps && (layers["hfs_wall"] < min_fraction_thin_wall * norm)
        layers["hfs_blanket"] -= (min_fraction_thin_wall * norm - layers["hfs_wall"])
        layers["hfs_wall"] = min_fraction_thin_wall * norm
    end
    if no_small_gaps && (layers["lfs_wall"] < min_fraction_thin_wall * norm)
        layers["lfs_blanket"] -= (min_fraction_thin_wall * norm - layers["lfs_wall"])
        layers["lfs_wall"] = min_fraction_thin_wall * norm
    end

    return layers
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
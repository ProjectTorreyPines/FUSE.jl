using DataFrames: DataFrames
using CSV: CSV

"""
    case_parameters(::Type{Val{:HDB5}}; tokamak::Union{String,Symbol}=:any, case=missing, database_case=missing)

For description of cases/variables see https://osf.io/593q6/
"""
function case_parameters(::Type{Val{:HDB5}}; tokamak::Union{String,Symbol}=:any, case=missing, database_case=missing)::Tuple{ParametersAllInits,ParametersAllActors}
    if !ismissing(database_case)
        data_row = load_hdb5(; database_case)
    elseif !ismissing(case)
        data_row = load_hdb5(tokamak)[case, :]
    else
        error("Specifcy either the case or database_case")
    end
    return case_parameters(data_row)
end

function case_parameters(data_row::DataFrames.DataFrameRow)
    ini = ParametersInits()
    act = ParametersActors()
    ini.general.casename = "HDB_$(data_row[:TOK])_$(data_row[:SHOT])"
    ini.general.init_from = :scalars

    # Equilibrium parameters
    ini.equilibrium.boundary_from = :scalars
    ini.equilibrium.B0 = data_row[:BT]
    ini.equilibrium.R0 = data_row[:RGEO]
    ini.equilibrium.Z0 = 0.0
    ini.equilibrium.ϵ = data_row[:AMIN] / data_row[:RGEO]
    ini.equilibrium.κ = data_row[:KAPPA]
    ini.equilibrium.δ = data_row[:DELTA]
    ini.equilibrium.ip = data_row[:IP]
    ini.equilibrium.pressure_core = IMAS.pressure_avg_from_beta_n(1.0, data_row[:AMIN], data_row[:BT], data_row[:IP]) * 3.0

    # Determine x-points
    if data_row[:CONFIG] == "SN"
        # upper single null
        ini.equilibrium.xpoints = :upper
    elseif data_row[:CONFIG] == "SN(L)"
        # lower single null
        ini.equilibrium.xpoints = :lower
    elseif data_row[:CONFIG] == "DN"
        # double null
        ini.equilibrium.xpoints = :double
    else
        # no x-points
        ini.equilibrium.xpoints = :none
    end

    # to match the experimental volume and area:
    mxhb = MXHboundary(ini; target_volume=data_row[:VOL], target_area=data_row[:AREA])
    mxh = IMAS.MXH(mxhb.r_boundary, mxhb.z_boundary, 4)
    ini.equilibrium.MXH_params = IMAS.flat_coeffs(mxh)
    ini.equilibrium.boundary_from = :MXH_params

    # Core_profiles parameters
    ini.core_profiles.ne_setting = :ne_line
    ini.core_profiles.ne_value = data_row[:NEL]
    ini.core_profiles.ne_shaping = 0.9
    ini.core_profiles.Te_shaping = 1.8
    ini.core_profiles.Ti_Te_ratio = 1.0
    ini.core_profiles.zeff = data_row[:ZEFF]
    ini.core_profiles.rot_core = 10e3
    ini.core_profiles.ngrid = 201
    ini.core_profiles.bulk = :D
    ini.core_profiles.impurity = :C

    # hcd
    if data_row[:PNBI] > 0
        resize!(ini.nb_unit, 1)
        ini.nb_unit[1].power_launched = data_row[:PNBI]
        if data_row[:ENBI] > 0
            ini.nb_unit[1].beam_energy = data_row[:ENBI]
        else
            ini.nb_unit[1].beam_energy = 100e3
        end
        ini.nb_unit[1].beam_mass = 2.0
        ini.nb_unit[1].toroidal_angle = 18.0 * deg # 18 degrees assumed like DIII-D
    end
    if data_row[:PECRH] > 0
        resize!(ini.ec_launcher, 1)
        ini.ec_launcher[1].power_launched = data_row[:PECRH]
    end
    if data_row[:PICRH] > 0
        resize!(ini.ic_antenna, 1)
        ini.ic_antenna[1].power_launched = data_row[:PICRH]
    end

    #### ACT ####

    act.ActorPedestal.density_match = :ne_line
    act.ActorFluxMatcher.evolve_pedestal = false

    act.ActorTGLF.tglfnn_model = "sat1_em_d3d"

    return ini, act
end

function load_hdb5(tokamak::Union{String,Symbol}=:all; maximum_ohmic_fraction::Float64=0.25, database_case::Union{Int,Missing}=missing, extra_signal_names=Union{String,Symbol}[])
    # For description of variables see https://osf.io/593q6/
    run_df = CSV.read(joinpath(__FUSE__, "sample", "HDB5_compressed.csv"), DataFrames.DataFrame)
    run_df[:, "database_case"] = collect(1:length(run_df[:, "TOK"]))

    if !ismissing(database_case)
        return run_df[run_df.database_case.==database_case, :]
    end

    signal_names = [
        "TOK", "SHOT", "TIME",
        "AMIN", "KAPPA", "DELTA",
        "NEL", "ZEFF", "TAUTH",
        "RGEO", "BT", "IP",
        "PNBI", "ENBI", "PICRH",
        "PECRH", "POHM", "MEFF",
        "VOL", "AREA", "WTH",
        "CONFIG"
    ]
    signal_names = vcat(signal_names, extra_signal_names)

    # subselect on the signals of interest
    run_df = run_df[:, signal_names]

    # only retain cases for which all signals have data
    run_df = run_df[DataFrames.completecases(run_df), :]

    # some basic filters
    run_df = run_df[(run_df.TOK.!="T10").&(run_df.TOK.!="TDEV").&(run_df.KAPPA.>1.0).&(run_df.DELTA.<0.79).&(1.6 .< run_df.MEFF .< 2.2).&(1.1 .< run_df.ZEFF .< 5.9), :]

    # Filter cases where the ohmic power is dominating
    run_df[:, "Paux"] = run_df[:, "PNBI"] .+ run_df[:, "PECRH"] .+ run_df[:, "PICRH"] .+ run_df[:, "POHM"]
    run_df = run_df[run_df[:, "POHM"].<maximum_ohmic_fraction.*(run_df[:, "Paux"].-run_df[:, "POHM"]), :]
    if Symbol(tokamak) ∉ (:all, :any)
        run_df = run_df[run_df.TOK.==String(tokamak), :]
    end

    return run_df
end

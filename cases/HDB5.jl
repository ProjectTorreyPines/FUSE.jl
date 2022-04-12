import DataFrames
import CSV

# For description of cases/variables see https://osf.io/593q6/
"""
    case_parameters(::Type{Val{:HDB5}}; tokamak::Union{String,Symbol}=:any, case::Integer,database_id::Integer)
"""
function case_parameters(::Type{Val{:HDB5}}; tokamak::Union{String,Symbol}=:any, case=missing, database_id=missing)
    if !ismissing(database_id)
        data_row = load_hdb5(database_id=database_id)
    elseif !ismissing(case)
        data_row = load_hdb5(tokamak)[case, :]
    else
        error("Specifcy either the case or database_id")
    end
    case_parameters(data_row)
end

function case_parameters(data_row::DataFrames.DataFrameRow)
    ini = InitParameters()
    act = ActorParameters()
    ini.general.casename = "HDB_$(data_row[:TOK])_$(data_row[:SHOT])"
    ini.general.init_from = :scalars

    # Equilibrium parameters
    ini.equilibrium.B0 = data_row[:BT]
    ini.equilibrium.R0 = data_row[:RGEO]
    ini.equilibrium.Z0 = 0.0
    ini.equilibrium.ϵ = data_row[:AMIN] / data_row[:RGEO]
    ini.equilibrium.κ = data_row[:KAPPA]
    ini.equilibrium.δ = data_row[:DELTA]
    ini.equilibrium.βn = 1.0
    ini.equilibrium.ip = data_row[:IP]

    act.SolovevActor.area = data_row[:AREA]
    act.SolovevActor.volume = data_row[:VOL]

    x_point = (data_row[:RGEO] * (1 - 1.1 * data_row[:DELTA] * data_row[:AMIN] / data_row[:RGEO]), data_row[:RGEO] * 1.1 * data_row[:KAPPA] * data_row[:AMIN] / data_row[:RGEO])

    # Determine x-points
    if data_row[:CONFIG] == "SN"
        # upper single null
        x_point = (x_point[1], x_point[2])
        symmetric = false
    elseif data_row[:CONFIG] == "SN(L)"
        # lower single null
        x_point = (x_point[1], -x_point[2])
        symmetric = false
    elseif data_row[:CONFIG] == "DN"
        # double null
        x_point = (x_point[1], x_point[2])
        symmetric = true
    else
        # no x-points
        x_point = false
        symmetric = true
    end

    ini.equilibrium.x_point = x_point
    ini.equilibrium.symmetric = symmetric

    # Core_profiles parameters
    ini.core_profiles.ne_ped = data_row[:NEL] / 1.3
    ini.core_profiles.n_peaking = 1.5
    ini.core_profiles.T_shaping = 1.8
    ini.core_profiles.w_ped = 0.03
    ini.core_profiles.zeff = data_row[:ZEFF]
    ini.core_profiles.rot_core = 50e3
    ini.core_profiles.ngrid = 201
    ini.core_profiles.bulk = :D
    ini.core_profiles.impurity = :C

    # nbi
    if data_row[:PNBI] > 0
        ini.nbi.power_launched = data_row[:PNBI]
        if data_row[:ENBI] > 0
            ini.nbi.beam_energy = data_row[:ENBI]
        else
            ini.nbi.beam_energy = 100e3
        end
        ini.nbi.beam_mass = 2
        ini.nbi.toroidal_angle = 0.0
    end

    if data_row[:PECRH] > 0
        ini.ec.power_launched = data_row[:PECRH]
    end

    if data_row[:PICRH] > 0
        ini.ic.power_launched = data_row[:PICRH]
    end

    return set_new_base!(ini), set_new_base!(act)
end

function load_hdb5(tokamak::T=:all; maximum_ohmic_fraction=0.25, database_id=missing, extra_signal_names=T[]) where {T<:Union{String,Symbol}}
    # For description of variables see https://osf.io/593q6/
    run_df = CSV.read(joinpath(dirname(abspath(@__FILE__)), "..", "sample", "HDB5_compressed.csv"), DataFrames.DataFrame)
    run_df[:,"database_id"] = collect(StepRange(1,1,length(run_df[:,"TOK"])))

    if !ismissing(database_id)
        return run_df[run_df.database_id .== database_id, :]
    end

    signal_names = ["TOK", "SHOT", "AMIN", "KAPPA", "DELTA", "NEL", "ZEFF", "TAUTH","AUXHEAT", "RGEO", "BT", "IP", "PNBI", "ENBI", "PICRH", "PECRH", "POHM", "MEFF", "VOL", "AREA", "WTH", "CONFIG"]
    signal_names = vcat(signal_names, extra_signal_names)
    # subselect on the signals of interest
    run_df = run_df[:, signal_names]
    # only retain cases for which all signals have data
    run_df = run_df[DataFrames.completecases(run_df), :]
    # some basic filters
    run_df = run_df[(run_df.TOK.!="T10").&(run_df.TOK.!="TDEV").&(run_df.KAPPA.>1.0).&(1.6 .< run_df.MEFF .< 2.2).&(1.1 .< run_df.ZEFF .< 5.9), :]
    # Filter cases where the ohmic power is dominating
    run_df[:,"Paux"] = run_df[:,"PNBI"] .+ run_df[:,"PECRH"] .+ run_df[:,"PICRH"] .+ run_df[:,"POHM"]
    run_df = run_df[run_df[:,"POHM"] .< maximum_ohmic_fraction .* (run_df[:,"Paux"] .- run_df[:,"POHM"]),:]
    if !(Symbol(tokamak) in [:all, :any])
        run_df = run_df[run_df.TOK.==String(tokamak), :]
    end
    return run_df
end
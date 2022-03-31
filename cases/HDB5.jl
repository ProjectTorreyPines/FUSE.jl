import DataFrames
import CSV

# For description of cases/variables see https://osf.io/593q6/
function Parameters(::Type{Val{:HDB5}}; tokamak::Union{String,Symbol}=:any, case::Integer)
    data_row = load_hdb5(tokamak)[case, :]
    Parameters(data_row)
end

function Parameters(data_row::DataFrames.DataFrameRow)
    par = Parameters()
    par.general.casename = "HDB_$(data_row[:TOK])_$(data_row[:SHOT]))"
    par.general.init_from = :scalars

    # Equilibrium parameters
    par.equilibrium.B0 = abs(data_row[:BT])
    par.equilibrium.R0 = data_row[:RGEO]
    par.equilibrium.Z0 = 0.0
    par.equilibrium.ϵ = data_row[:AMIN] / data_row[:RGEO]
    par.equilibrium.κ = data_row[:KAPPA]
    par.equilibrium.δ = data_row[:DELTA]
    par.equilibrium.βn = 1.0
    par.equilibrium.area = data_row[:AREA]
    par.equilibrium.volume = data_row[:VOL]
    par.equilibrium.ip = abs(data_row[:IP])
    par.equilibrium.x_point = false
    par.equilibrium.symmetric = true

    # Core_profiles parameters
    par.core_profiles.ne_ped = data_row[:NEL] / 1.3
    par.core_profiles.n_peaking = 1.5
    par.core_profiles.T_shaping = 1.8
    par.core_profiles.w_ped = 0.03
    par.core_profiles.zeff = data_row[:ZEFF]
    par.core_profiles.rot_core = 50e3
    par.core_profiles.ngrid = 201
    par.core_profiles.bulk = :D
    par.core_profiles.impurity = :C

    # nbi
    if data_row[:PNBI] > 0
        par.nbi.power_launched = data_row[:PNBI]
        if data_row[:ENBI] > 0
            par.nbi.beam_energy = data_row[:ENBI]
        else
            par.nbi.beam_energy = 100e3
        end
        par.nbi.beam_mass = 2
        par.nbi.toroidal_angle = 0.0
    end

    # ohmic 

    if data_row[:POHM] > 0
        par.oh.ohmic_heating = data_row[:POHM]
    end

    if data_row[:PECRH] > 0
        par.ec.power_launched = data_row[:PECRH]
    end

    if data_row[:PICRH] > 0
        par.ic.power_launched = data_row[:PICRH]
    end

    return par
end

function load_hdb5(tokamak::T=:all, extra_signal_names=T[]) where {T<:Union{String,Symbol}}
    # Set up the database to run
    # For description of variables see https://osf.io/593q6/
    run_df = CSV.read(joinpath(dirname(abspath(@__FILE__)), "..", "sample", "HDB5_compressed.csv"), DataFrames.DataFrame)
    signal_names = ["TOK", "SHOT", "AMIN", "KAPPA", "DELTA", "NEL", "ZEFF", "TAUTH", "RGEO", "BT", "IP", "PNBI", "ENBI", "PICRH", "PECRH", "POHM", "MEFF", "VOL", "AREA", "WTH", "CONFIG"]
    signal_names = vcat(signal_names, extra_signal_names)
    # subselect on the signals of interest
    run_df = run_df[:, signal_names]
    # only retain cases for which all signals have data
    run_df = run_df[DataFrames.completecases(run_df), :]
    # some basic filters
    run_df = run_df[(run_df.TOK.!="T10").&(run_df.TOK.!="TDEV").&(run_df.KAPPA.>1.0).&(1.6 .< run_df.MEFF .< 2.2).&(1.1 .< run_df.ZEFF .< 5.9), :]
    if !(Symbol(tokamak) in [:all,:any])
        run_df = run_df[run_df.TOK.==String(tokamak), :]
    end
    return run_df
end                                                                                                                                                                                                   
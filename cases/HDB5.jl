import DataFrames
import CSV

function Parameters(::Type{Val{:HDB5}}, filename::String, case_number::Integer)
    run_df = CSV.read(filename, DataFrames.DataFrame)
    data_row = run_df[case_number, :]
    Parameters(:HDB5; data_row)
end

function Parameters(::Type{Val{:HDB5}}; data_row::Union{NamedTuple,AbstractDict})
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
    if data_row[:PNBI] + data_row[:POHM] > 0.0
        par.nbi.beam_power = data_row[:PNBI] + data_row[:POHM]
        if data_row[:PNBI] == 0.0
            par.nbi.beam_energy = 1e9
        else
            par.nbi.beam_energy = data_row[:ENBI]
        end
        par.nbi.beam_mass = 2
        par.nbi.toroidal_angle = 0.0
    end

    if data_row[:PECRH] > 0
        par.ec.power_launched = data_row[:PECRH]
    end

    if data_row[:PICRH] > 0
        par.ic.power_launched = data_row[:PICRH]
    end

    return par
end
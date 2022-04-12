import DataFrames
import CSV

# For description of cases/variables see https://osf.io/593q6/
```
    Parameters(::Type{Val{:HDB5}}; tokamak::Union{String,Symbol}=:any, case::Integer,database_id::Integer)
```
function Parameters(::Type{Val{:HDB5}}; tokamak::Union{String,Symbol}=:any; case::Integer=missing, database_id::Integer=missing)
    if !ismissing(database_id)
        data_row = load_hdb5(database_id=database_id)
    elseif !ismissing(case)
        data_row = load_hdb5(tokamak)[case_prefilter, :]
    else
        error("Specifcy either the case or case_prefilter")
    end
    Parameters(data_row)
end

function Parameters(data_row::DataFrames.DataFrameRow)
    par = Parameters()
    par.general.casename = "HDB_$(data_row[:TOK])_$(data_row[:SHOT])"
    par.general.init_from = :scalars

    # Equilibrium parameters
    par.equilibrium.B0 = data_row[:BT]
    par.equilibrium.R0 = data_row[:RGEO]
    par.equilibrium.Z0 = 0.0
    par.equilibrium.ϵ = data_row[:AMIN] / data_row[:RGEO]
    par.equilibrium.κ = data_row[:KAPPA]
    par.equilibrium.δ = data_row[:DELTA]
    par.equilibrium.βn = 1.0
    par.equilibrium.area = data_row[:AREA]
    par.equilibrium.volume = data_row[:VOL]
    par.equilibrium.ip = data_row[:IP]

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

    par.equilibrium.x_point = x_point
    par.equilibrium.symmetric = symmetric

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

    if data_row[:PECRH] > 0
        par.ec.power_launched = data_row[:PECRH]
    end

    if data_row[:PICRH] > 0
        par.ic.power_launched = data_row[:PICRH]
    end
    show(data_row,allcols=true)
    return par
end

function load_hdb5(tokamak::T=:all, extra_signal_names=T[]) where {T<:Union{String,Symbol};maximum_ohmic_fraction=0.25; database_id::Integer=missing}
    # Set up the database to run
    # For description of variables see https://osf.io/593q6/
    run_df = CSV.read(joinpath(dirname(abspath(@__FILE__)), "..", "sample", "HDB5_compressed.csv"), DataFrames.DataFrame)
    run_df[:,"database_id"] = collect(StepRange(0,1,len(run_df[:,"TOK"])))
    if !missing(case_prefilter)
        return run_df[run_df.database_id = database_id, :]
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


#pth08 = 0.049 * (abs(bt) ** 0.8) * (ne ** 0.72) * (surf_area ** 0.94)  # MW
# P_TH in MW for ne in m^20 m^-3, Bt in T, and area in m^2

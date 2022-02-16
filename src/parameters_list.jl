struct ParametersKwargs
    group::Symbol
    kw
    not_used
end

function ParametersKwargs(group, kw)
    ParametersKwargs(group, Dict(key => value for (key, value) in kw), Dict(key => value for (key, value) in kw))
end

function Base.getindex(pa::ParametersKwargs, key)
    if key in keys(pa.kw)
        if key in keys(pa.not_used)
            pop!(pa.not_used, key)
        end
        return pa.kw[key]
    else
        error("Need to set keyword argument `$key` when calling Parameters($(repr(pa.group)); ...)")
    end
end

function Base.pop!(pa::ParametersKwargs, key)
    if key in keys(pa.not_used)
        pop!(pa.not_used, key)
    end
    return pop!(pa.kw, key)
end

function Base.length(pa::ParametersKwargs)
    return length(pa.kw)
end

top_level_parameters = [:general, :equilibrium, :core_profiles, :coil, :build, :gasc, :ods]

function Parameters()
    params = Parameters(Dict{Symbol,Union{Parameter,Parameters}}())
    for item in top_level_parameters
        setproperty!(params, item, Parameters(item))
    end
    return params
end

function Parameters(group::Symbol; kw...)
    kw = ParametersKwargs(group, kw)
    if group in top_level_parameters
        params = Parameters(Dict{Symbol,Union{Parameter,Parameters}}())

        if group == :general
            general = params
            options = [
                :ods => "Load data from ODS saved in .json format",
                :scalars => "Initialize FUSE run form scalar FUSE parameters",
                :gasc => "Initialize FUSE run form GASC output file saved in .json format"]
            general.init_from = Switch(options, "", "Initialize run from")

        elseif group == :equilibrium
            equilibrium = params
            equilibrium.B0 = Entry(Real, IMAS.equilibrium__vacuum_toroidal_field, :b0)
            equilibrium.R0 = Entry(Real, IMAS.equilibrium__vacuum_toroidal_field, :r0)
            equilibrium.Z0 = Entry(Real, "m", "Z offset of the machine midplane"; default = 0.0)
            equilibrium.ϵ = Entry(Real, "", "Plasma aspect ratio")
            equilibrium.δ = Entry(Real, IMAS.equilibrium__time_slice___boundary, :triangularity)
            equilibrium.κ = Entry(Real, IMAS.equilibrium__time_slice___boundary, :elongation)
            equilibrium.βn = Entry(Real, IMAS.equilibrium__time_slice___global_quantities, :beta_normal)
            equilibrium.ip = Entry(Real, IMAS.equilibrium__time_slice___global_quantities, :ip)
            equilibrium.x_point = Entry(Bool, IMAS.equilibrium__time_slice___boundary, :x_point)
            equilibrium.symmetric = Entry(Bool, "", "Is plasma up-down symmetric")
            equilibrium.ngrid = Entry(Int, "", "Resolution of the equilibrium grid"; default = 129)
            equilibrium.field_null_surface = Entry(Real, "", "ψn value of the field_null_surface. Disable with 0.0"; default = 0.25)#, min=0.0, max=1.0)

        elseif group == :core_profiles
            core_profiles = params
            core_profiles.ne_ped = Entry(Real, "m^-3", "Pedestal electron density")
            core_profiles.ne_peaking = Entry(Real, "", "Electron density peaking")
            core_profiles.Te_ped = Entry(Real, "eV", "Pedestal electron temperature")
            core_profiles.Te_peaking = Entry(Real, "", "Electron temperature peaking")
            core_profiles.w_ped = Entry(Real, "", "Pedestal width expressed in fraction of ψₙ")
            core_profiles.zeff = Entry(Real, "", "Effective ion charge")
            core_profiles.P_co_nbi = Entry(Real, "W", "Co-current NBI input power")
            core_profiles.ngrid = Entry(Int, "", "Resolution of the core_profiles grid"; default = 101)
            core_profiles.bulk = Entry(Symbol, "", "Bulk ion species")
            core_profiles.impurity = Entry(Symbol, "", "Impurity ion species")

        elseif group == :coil
            coil = params
            options = [
                :point => "one filament per coil",
                :simple => "like :point, but OH coils have three filaments",
                :corners => "like :simple, but PF coils have filaments at the four corners",
                :realistic => "hundreds of filaments per coil (very slow!)"]
            coil.green_model = Switch(options, "", "Model to be used for the Greens function table of the PF coils"; default = :simple)

        elseif group == :build
            build = params
            build.n_oh_coils = Entry(Int, "", "Number of OH coils")
            build.n_pf_coils_inside = Entry(Int, "", "Number of PF coils inside of the TF")
            build.n_pf_coils_outside = Entry(Int, "", "Number of PF coils outside of the TF")

            build.is_nuclear_facility = Entry(Bool, "", "Is this a nuclear facility")

        elseif group == :gasc
            gasc = params
            gasc.filename = Entry(String, "", "Output GASC .json file from which data will be loaded")
            gasc.case = Entry(Int, "", "Number of the GASC run to load")
            gasc.no_small_gaps = Entry(Bool, "", "Remove small gaps from the GASC radial build"; default = true)

        elseif group == :ods
            ods = params
            ods.filename = Entry(String, "", "ODS.json file from which equilibrium is loaded")
        end

    else
        params = Parameters()

        if group == :ITER
            params.general.init_from = kw[:init_from]

            if kw[:init_from] == :ods
                params.ods.filename = joinpath(dirname(abspath(@__FILE__)), "..", "sample", "ITER_eq_ods.json")
            else
                params.equilibrium.R0 = 6.2
                params.equilibrium.ϵ = 0.32
                params.equilibrium.κ = 1.85
                params.equilibrium.δ = 0.485
                params.equilibrium.B0 = 5.3
                params.equilibrium.Z0 = 0.4
                params.equilibrium.ip = 15.E6
                params.equilibrium.βn = 2.0
                params.equilibrium.x_point = true
                params.equilibrium.symmetric = false
            end

            params.build.is_nuclear_facility = true
            params.build.n_oh_coils = 6
            params.build.n_pf_coils_inside = 0
            params.build.n_pf_coils_outside = 6

            params.core_profiles.ne_ped = 7E19
            params.core_profiles.ne_peaking = 1.5
            params.core_profiles.Te_ped = 1E3
            params.core_profiles.Te_peaking = 3
            params.core_profiles.w_ped = 0.08
            params.core_profiles.zeff = 2.5
            params.core_profiles.P_co_nbi = 20e6
            params.core_profiles.bulk = :DT
            params.core_profiles.impurity = :Ne

        elseif group == :CAT
            params.general.init_from = :ods

            params.ods.filename = joinpath(dirname(abspath(@__FILE__)), "..", "sample", "CAT_eq_ods.json")

            params.build.is_nuclear_facility = false
            params.build.n_oh_coils = 6
            params.build.n_pf_coils_inside = 0
            params.build.n_pf_coils_outside = 6

            params.core_profiles.ne_ped = 7E19
            params.core_profiles.ne_peaking = 1.5
            params.core_profiles.Te_ped = 1E3
            params.core_profiles.Te_peaking = 3
            params.core_profiles.w_ped = 0.08
            params.core_profiles.zeff = 2.5
            params.core_profiles.P_co_nbi = 20e6
            params.core_profiles.bulk = :DT
            params.core_profiles.impurity = :Ne

        elseif group == :D3D
            params.general.init_from = :ods

            params.ods.filename = joinpath(dirname(abspath(@__FILE__)), "..", "sample", "D3D_eq_ods.json")

            params.build.is_nuclear_facility = false
            params.build.n_oh_coils = 20
            params.build.n_pf_coils_inside = 18
            params.build.n_pf_coils_outside = 0

            params.core_profiles.ne_ped = 5E19
            params.core_profiles.ne_peaking = 1.5
            params.core_profiles.Te_ped = 500
            params.core_profiles.Te_peaking = 3
            params.core_profiles.w_ped = 0.08
            params.core_profiles.zeff = 2.0
            params.core_profiles.P_co_nbi = 5e6
            params.core_profiles.bulk = :D
            params.core_profiles.impurity = :C

        elseif group == :FPP
            params.general.init_from = :gasc

            params.gasc.filename = joinpath(dirname(abspath(@__FILE__)), "..", "sample", "FPP_fBS_PBpR_scan.json")
            params.gasc.case = 59
            params.gasc.no_small_gaps = true

            params.build.is_nuclear_facility = true
            params.build.n_oh_coils = 6
            params.build.n_pf_coils_inside = 0
            params.build.n_pf_coils_outside = 6

            params.core_profiles.bulk = :DT
            params.core_profiles.impurity = :Ne
        else
            throw(InexistentParameterException(group))
        end
        if length(kw.not_used) > 0
            error("Parameters($(repr(group))) did not use the following arguments parameters: $(collect(keys(kw.kw)))")
        end
        set_new_base(params)
    end

    return params
end

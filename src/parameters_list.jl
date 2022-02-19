function Parameters()
    par = Parameters(Dict{Symbol,Union{Parameter,Parameters}}())
    for item in [:general, :equilibrium, :core_profiles, :pf_active, :nbi, :build, :gasc, :ods]
        setproperty!(par, item, Parameters(item))
    end
    return par
end

parametersDispatcher = Dict()

function Parameters(group::Symbol; kw...)
    if group in keys(parametersDispatcher)
        return Parameters(parametersDispatcher[group]; kw...)
    end

    par = Parameters(Dict{Symbol,Union{Parameter,Parameters}}())

    if group == :general
        general = par
        options = [
            :ods => "Load data from ODS saved in .json format",
            :scalars => "Initialize FUSE run form scalar FUSE parameters",
            :gasc => "Initialize FUSE run form GASC output file saved in .json format"]
        general.init_from = Switch(options, "", "Initialize run from")

    elseif group == :equilibrium
        equilibrium = par
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
        core_profiles = par
        core_profiles.ne_ped = Entry(Real, "m^-3", "Pedestal electron density")
        core_profiles.n_peaking = Entry(Real, "", "Ratio of core/pedestal densities")
        core_profiles.T_shaping = Entry(Real, "", "Temperature shaping factor")
        core_profiles.w_ped = Entry(Real, "", "Pedestal width expressed in fraction of ψₙ")
        core_profiles.zeff = Entry(Real, "", "Effective ion charge")
        core_profiles.rot_core = Entry(Real, IMAS.core_profiles__profiles_1d, :rotation_frequency_tor_sonic)
        core_profiles.ngrid = Entry(Int, "", "Resolution of the core_profiles grid"; default = 101)
        core_profiles.bulk = Entry(Symbol, "", "Bulk ion species")
        core_profiles.impurity = Entry(Symbol, "", "Impurity ion species")

    elseif group == :pf_active
        pf_active = par
        options = [
            :point => "one filament per coil",
            :simple => "like :point, but OH coils have three filaments",
            :corners => "like :simple, but PF coils have filaments at the four corners",
            :realistic => "hundreds of filaments per coil (very slow!)"]
        pf_active.green_model = Switch(options, "", "Model to be used for the Greens function table of the PF coils"; default = :simple)
        pf_active.n_oh_coils = Entry(Int, "", "Number of OH coils")
        pf_active.n_pf_coils_inside = Entry(Int, "", "Number of PF coils inside of the TF")
        pf_active.n_pf_coils_outside = Entry(Int, "", "Number of PF coils outside of the TF")

    elseif group == :nbi
        nbi = par
        nbi.beam_power = Entry(Union{Real,Vector{Real}}, "W", "Beam power")
        nbi.beam_energy = Entry(Union{Real,Vector{Real}}, "eV", "Beam energy")
        nbi.beam_mass = Entry(Union{Real,Vector{Real}}, "AU", "Beam mass"; default = 2.0)
        nbi.toroidal_angle = Entry(Union{Real,Vector{Real}}, "rad", "toroidal angle of injection"; default = 0.0)

    elseif group == :build
        build = par
        build.is_nuclear_facility = Entry(Bool, "", "Is this a nuclear facility")

    elseif group == :gasc
        gasc = par
        gasc.filename = Entry(String, "", "Output GASC .json file from which data will be loaded")
        gasc.case = Entry(Int, "", "Number of the GASC run to load")
        gasc.no_small_gaps = Entry(Bool, "", "Remove small gaps from the GASC radial build"; default = true)

    elseif group == :ods
        ods = par
        ods.filename = Entry(String, "", "ODS.json file from which equilibrium is loaded")
    else
        throw(InexistentParameterException(group))
    end

    return par
end

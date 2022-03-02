using FusionMaterials: FusionMaterials

case_parameters = Symbol[]
for filename in readdir(joinpath(dirname(@__FILE__), "..", "cases"))
    push!(case_parameters, Symbol(splitext(filename)[1]))
    include("../cases/" * filename)
end

function Parameters()
    par = Parameters(Symbol[], Dict{Symbol,Union{Parameter,Parameters}}())
    for item in [:general, :equilibrium, :core_profiles, :pf_active, :oh, :tf, :nbi, :ec, :ic, :lh, :build, :gasc, :ods, :material]
        setproperty!(par, item, Parameters(item))
    end
    return par
end

function Parameters(::Type{Val{:general}})
    general = Parameters(nothing)
    general.casename = Entry(String, "", "Mnemonic name of the case being run")
    options = [
        :ods => "Load data from ODS saved in .json format",
        :scalars => "Initialize FUSE run form scalar FUSE parameters",
        :gasc => "Initialize FUSE run form GASC output file saved in .json format",
    ]
    general.init_from = Switch(options, "", "Initialize run from")
    return general
end

function Parameters(::Type{Val{:material}})
    material = Parameters(nothing)
    material.wall = Switch(FusionMaterials.available_materials("wall_materials"), "", "Material used for the wall"; default = "Steel, Stainless 316")
    material.blanket = Switch(FusionMaterials.available_materials("blanket_materials"), "", "Material used for blanket coils")
    material.shield = Switch(FusionMaterials.available_materials("shield_materials"), "", "Material used for the shield")
    return material
end

function Parameters(::Type{Val{:equilibrium}})
    equilibrium = Parameters(nothing)
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
    return equilibrium
end

function Parameters(::Type{Val{:core_profiles}})
    core_profiles = Parameters(nothing)
    core_profiles.ne_ped = Entry(Real, "m^-3", "Pedestal electron density")
    core_profiles.n_peaking = Entry(Real, "", "Ratio of core/pedestal densities")
    core_profiles.T_shaping = Entry(Real, "", "Temperature shaping factor")
    core_profiles.w_ped = Entry(Real, "", "Pedestal width expressed in fraction of ψₙ")
    core_profiles.zeff = Entry(Real, "", "Effective ion charge")
    core_profiles.rot_core = Entry(Real, IMAS.core_profiles__profiles_1d, :rotation_frequency_tor_sonic)
    core_profiles.ngrid = Entry(Int, "", "Resolution of the core_profiles grid"; default = 101)
    core_profiles.bulk = Entry(Symbol, "", "Bulk ion species")
    core_profiles.impurity = Entry(Symbol, "", "Impurity ion species")
    return core_profiles
end

function Parameters(::Type{Val{:pf_active}})
    pf_active = Parameters(nothing)
    options = [
        :point => "one filament per coil",
        :simple => "like :point, but OH coils have three filaments",
        :corners => "like :simple, but PF coils have filaments at the four corners",
        :realistic => "hundreds of filaments per coil (very slow!)",
    ]
    pf_active.green_model = Switch(options, "", "Model used for the Greens function calculation"; default = :simple)
    pf_active.n_oh_coils = Entry(Int, "", "Number of OH coils")
    pf_active.n_pf_coils_inside = Entry(Int, "", "Number of PF coils inside of the TF")
    pf_active.n_pf_coils_outside = Entry(Int, "", "Number of PF coils outside of the TF")
    pf_active.technology = Entry(Parameters, "", "PF coil technology")
    return pf_active
end

function Parameters(::Type{Val{:tf}})
    tf = Parameters(nothing)
    tf.n_coils = Entry(Int, "", "Number of TF coils")
    options = [1 => "PricetonD", 2 => "Rectangle", 3 => "TrippleArc", 4 => "Miller", 5 => "Spline"]
    #        tf.shape = Switch(options, "", "Shape of the TF coils"; default=:TrippleArc)
    tf.shape = Entry(Int, "", "Shape of the TF coils"; default = 3)
    tf.technology = Entry(Parameters, "", "TF coil technology")
    return tf
end

function Parameters(::Type{Val{:oh}})
    oh = Parameters(nothing)
    oh.technology = Entry(Parameters, "", "OH coil technology")
    return oh
end

function Parameters(::Type{Val{:nbi}})
    nbi = Parameters(nothing)
    nbi.beam_power = Entry(Union{X,Vector{X}} where {X<:Real}, "W", "Beam power")
    nbi.beam_energy = Entry(Union{X,Vector{X}} where {X<:Real}, "eV", "Beam energy")
    nbi.beam_mass = Entry(Union{X,Vector{X}} where {X<:Real}, "AU", "Beam mass"; default = 2.0)
    nbi.toroidal_angle = Entry(Union{X,Vector{X}} where {X<:Real}, "rad", "toroidal angle of injection"; default = 0.0)
    return nbi
end

function Parameters(::Type{Val{:ec}})
    ec = Parameters(nothing)
    ec.power_launched = Entry(Union{X,Vector{X}} where {X<:Real}, "W", "EC launched power")
    return ec
end

function Parameters(::Type{Val{:ic}})
    ic = Parameters(nothing)
    ic.power_launched = Entry(Union{X,Vector{X}} where {X<:Real}, "W", "IC launched power")
    return ic
end

function Parameters(::Type{Val{:lh}})
    lh = Parameters(nothing)
    lh.power_launched = Entry(Union{X,Vector{X}} where {X<:Real}, "W", "LH launched power")
    return lh
end

function Parameters(::Type{Val{:build}})
    build = Parameters(nothing)
    build.is_nuclear_facility = Entry(Bool, "", "Is this a nuclear facility")
    build.symmetric = Entry(Bool, "", "Is the build up-down symmetric")
    return build
end

function Parameters(::Type{Val{:gasc}})
    gasc = Parameters(nothing)
    gasc.filename = Entry(String, "", "Output GASC .json file from which data will be loaded")
    gasc.case = Entry(Int, "", "Number of the GASC run to load")
    gasc.no_small_gaps = Entry(Bool, "", "Remove small gaps from the GASC radial build"; default = true)
    return gasc
end

function Parameters(::Type{Val{:ods}})
    ods = Parameters(nothing)
    ods.filename = Entry(String, "", "ODS.json file from which equilibrium is loaded")
    return ods
end

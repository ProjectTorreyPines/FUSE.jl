using FusionMaterials: FusionMaterials

"""
    ParametersInit()

Generates initalization parameters 
"""
function ParametersInit()
    par = ParametersInit(missing, WeakRef(missing), Dict{Symbol,Union{Parameter,ParametersInit}}())
    for item in [:general, :equilibrium, :core_profiles, :pf_active, :oh, :tf, :center_stack, :nbi, :ec_launchers, :ic_antennas, :lh_antennas, :build, :gasc, :ods, :material]
        setproperty!(par, item, ParametersInit(item))
    end
    return par
end

function ParametersInit(::Type{Val{:general}})
    general = ParametersInit(nothing)
    general.casename = Entry(String, "", "Mnemonic name of the case being run")
    options = [
        :ods => "Load data from ODS saved in .json format (where possible, and fallback on scalars otherwise)",
        :scalars => "Initialize FUSE run from scalar parameters"
    ]
    general.init_from = Switch(options, "", "Initialize run from")
    return general
end

function ParametersInit(::Type{Val{:material}})
    material = ParametersInit(nothing)
    material.wall = Switch(FusionMaterials.available_materials("wall_materials"), "", "Material used for the wall"; default="Steel, Stainless 316")
    material.blanket = Switch(FusionMaterials.available_materials("blanket_materials"), "", "Material used for blanket coils")
    material.shield = Switch(FusionMaterials.available_materials("shield_materials"), "", "Material used for the shield")
    return material
end

function ParametersInit(::Type{Val{:equilibrium}})
    equilibrium = ParametersInit(nothing)
    equilibrium.B0 = Entry(Real, IMAS.equilibrium__vacuum_toroidal_field, :b0)
    equilibrium.R0 = Entry(Real, IMAS.equilibrium__vacuum_toroidal_field, :r0)
    equilibrium.Z0 = Entry(Real, "m", "Z offset of the machine midplane"; default=0.0)
    equilibrium.ϵ = Entry(Real, "", "Plasma aspect ratio")
    equilibrium.δ = Entry(Real, IMAS.equilibrium__time_slice___boundary, :triangularity)
    equilibrium.κ = Entry(Real, IMAS.equilibrium__time_slice___boundary, :elongation)
    equilibrium.βn = Entry(Real, IMAS.equilibrium__time_slice___global_quantities, :beta_normal)
    equilibrium.ip = Entry(Real, IMAS.equilibrium__time_slice___global_quantities, :ip)
    equilibrium.x_point = Entry(Union{NTuple{2},Bool}, IMAS.equilibrium__time_slice___boundary, :x_point)
    equilibrium.symmetric = Entry(Bool, "", "Is plasma up-down symmetric")
    equilibrium.ngrid = Entry(Int, "", "Resolution of the equilibrium grid"; default=129)
    equilibrium.field_null_surface = Entry(Real, "", "ψn value of the field_null_surface. Disable with 0.0"; default=0.25)#, min=0.0, max=1.0)
    return equilibrium
end

function ParametersInit(::Type{Val{:core_profiles}})
    core_profiles = ParametersInit(nothing)
    core_profiles.ne_ped = Entry(Real, "m^-3", "Pedestal electron density")
    core_profiles.greenwald_fraction = Entry(Real, "", "Greenwald fraction, ne_vol / ne_gw")
    core_profiles.T_shaping = Entry(Real, "", "Temperature shaping factor")
    core_profiles.w_ped = Entry(Real, "", "Pedestal width expressed in fraction of ψₙ")
    core_profiles.zeff = Entry(Real, "", "Effective ion charge")
    core_profiles.rot_core = Entry(Real, IMAS.core_profiles__profiles_1d, :rotation_frequency_tor_sonic)
    core_profiles.ngrid = Entry(Int, "", "Resolution of the core_profiles grid"; default=101)
    core_profiles.bulk = Entry(Symbol, "", "Bulk ion species")
    core_profiles.impurity = Entry(Symbol, "", "Impurity ion species")
    core_profiles.helium_fraction = Entry(Real, "", "Helium density / electron density fraction")
    core_profiles.ejima = Entry(Real, "", "Ejima coefficient"; default=0.4)
    return core_profiles
end

function ParametersInit(::Type{Val{:pf_active}})
    pf_active = ParametersInit(nothing)
    pf_active.n_oh_coils = Entry(Int, "", "Number of OH coils")
    pf_active.n_pf_coils_inside = Entry(Int, "", "Number of PF coils inside of the TF")
    pf_active.n_pf_coils_outside = Entry(Int, "", "Number of PF coils outside of the TF")
    pf_active.technology = ParametersInit(:coil_technology)
    return pf_active
end

function ParametersInit(::Type{Val{:tf}})
    tf = ParametersInit(nothing)
    tf.n_coils = Entry(Int, "", "Number of TF coils")
    options = [:princeton_D_exact, :princeton_D, :princeton_D_scaled, :rectangle, :triple_arc, :miller, :spline]
    tf.shape = Switch(options, "", "Shape of the TF coils"; default=:princeton_D_scaled)
    tf.ripple = Entry(Real, "", "Fraction of toroidal field ripple evaluated at the outermost radius of the plasma chamber"; default=0.01)
    tf.technology = ParametersInit(:coil_technology)
    return tf
end

function ParametersInit(::Type{Val{:oh}})
    oh = ParametersInit(nothing)
    oh.technology = ParametersInit(:coil_technology)
    oh.flattop_duration = Entry(Real, "s", "Duration of the flattop (use Inf for steady-state)")
    return oh
end

function ParametersInit(::Type{Val{:center_stack}})
    center_stack = ParametersInit(nothing)
    center_stack.bucked = Entry(Bool, "", "flag for bucked boundary conditions between TF and OH (and center plug, if present"; default=false)
    center_stack.noslip = Entry(Bool, "", "flag for no slip conditions between TF and OH (and center plug, if present)"; default=false)
    center_stack.plug = Entry(Bool, "", "flag for center plug"; default=false)
    return center_stack
end


function ParametersInit(::Type{Val{:nbi}})
    nbi = ParametersInit(nothing)
    nbi.power_launched = Entry(Union{X,Vector{X}} where {X<:Real}, "W", "Beam power")
    nbi.beam_energy = Entry(Union{X,Vector{X}} where {X<:Real}, "eV", "Beam energy")
    nbi.beam_mass = Entry(Union{X,Vector{X}} where {X<:Real}, "AU", "Beam mass"; default=2.0)
    nbi.toroidal_angle = Entry(Union{X,Vector{X}} where {X<:Real}, "rad", "toroidal angle of injection"; default=0.0)
    nbi.efficiency_conversion = Entry(Union{X,Vector{X}} where {X<:Real}, IMAS.nbi__unit___efficiency, :conversion)
    nbi.efficiency_transmission = Entry(Union{X,Vector{X}} where {X<:Real}, IMAS.nbi__unit___efficiency, :transmission)
    return nbi
end

function ParametersInit(::Type{Val{:ec_launchers}})
    ec_launchers = ParametersInit(nothing)
    ec_launchers.power_launched = Entry(Union{X,Vector{X}} where {X<:Real}, "W", "EC launched power")
    ec_launchers.efficiency_conversion = Entry(Union{X,Vector{X}} where {X<:Real}, IMAS.ec_launchers__launcher___efficiency, :conversion)
    ec_launchers.efficiency_transmission = Entry(Union{X,Vector{X}} where {X<:Real}, IMAS.ec_launchers__launcher___efficiency, :transmission)
    return ec_launchers
end

function ParametersInit(::Type{Val{:ic_antennas}})
    ic_antennas = ParametersInit(nothing)
    ic_antennas.power_launched = Entry(Union{X,Vector{X}} where {X<:Real}, "W", "IC launched power")
    ic_antennas.efficiency_conversion = Entry(Union{X,Vector{X}} where {X<:Real}, IMAS.ic_antennas__antenna___efficiency, :conversion)
    ic_antennas.efficiency_transmission = Entry(Union{X,Vector{X}} where {X<:Real}, IMAS.ic_antennas__antenna___efficiency, :transmission)
    ic_antennas.efficiency_coupling = Entry(Union{X,Vector{X}} where {X<:Real}, IMAS.ic_antennas__antenna___efficiency, :coupling)
    return ic_antennas
end

function ParametersInit(::Type{Val{:lh_antennas}})
    lh_antennas = ParametersInit(nothing)
    lh_antennas.power_launched = Entry(Union{X,Vector{X}} where {X<:Real}, "W", "LH launched power")
    lh_antennas.efficiency_conversion = Entry(Union{X,Vector{X}} where {X<:Real}, IMAS.lh_antennas__antenna___efficiency, :conversion)
    lh_antennas.efficiency_transmission = Entry(Union{X,Vector{X}} where {X<:Real}, IMAS.lh_antennas__antenna___efficiency, :transmission)
    lh_antennas.efficiency_coupling = Entry(Union{X,Vector{X}} where {X<:Real}, IMAS.lh_antennas__antenna___efficiency, :coupling)
    return lh_antennas
end

function ParametersInit(::Type{Val{:build}})
    build = ParametersInit(nothing)
    build.layers = Entry(DataStructures.OrderedDict, "m", "Sorted dictionary of layers thicknesses in radial build")
    build.blanket = Entry(Float64, "", "Fraction of blanket in radial build")
    build.shield = Entry(Float64, "", "Fraction of shield in radial build")
    build.vessel = Entry(Float64, "", "Fraction of vessel in radial build")
    build.symmetric = Entry(Bool, "", "Is the build up-down symmetric")
    build.n_first_wall_conformal_layers = Entry(Integer, "", "Number of layers that are conformal to the first wall"; default=1)
    return build
end

function ParametersInit(::Type{Val{:gasc}})
    gasc = ParametersInit(nothing)
    gasc.filename = Entry(String, "", "Output GASC .json file from which data will be loaded")
    gasc.case = Entry(Int, "", "Number of the GASC run to load")
    return gasc
end

function ParametersInit(::Type{Val{:ods}})
    ods = ParametersInit(nothing)
    ods.filename = Entry(String, "", "ODS.json file from which equilibrium is loaded")
    return ods
end

function ParametersInit(::Type{Val{:coil_technology}})
    coil_tech = ParametersInit(nothing)
    coil_tech.material = Switch(FusionMaterials.available_materials("magnet_materials"), "", "Technology used for the coil.")
    coil_tech.temperature = Entry(Real, "K", "Coil temperature")
    coil_tech.thermal_strain = Entry(Real, "", "Fraction of thermal expansion strain over maximum total strain on coil")
    coil_tech.JxB_strain = Entry(Real, "", "Fraction of maximum JxB strain over maximum total strain on coil")
    coil_tech.fraction_stainless = Entry(Real, "", "Fraction of stainless steel in the coil cross-sectional areas")
    coil_tech.fraction_void = Entry(Real, "", "Fraction of `void` in the coil cross-sectional area. Void is everything (like coolant) that is not structural nor conductor.")
    coil_tech.ratio_SC_to_copper = Entry(Real, "", "Fraction of superconductor to copper cross-sectional areas")
    return coil_tech
end

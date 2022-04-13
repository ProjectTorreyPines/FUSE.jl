using FusionMaterials: FusionMaterials

"""
    InitParameters()

Generates initalization parameters 
"""
function InitParameters()
    par = InitParameters(Symbol[], Dict{Symbol,Union{Parameter,InitParameters}}())
    for item in [:general, :equilibrium, :core_profiles, :pf_active, :oh, :tf, :center_stack, :nbi, :ec, :ic, :lh, :build, :gasc, :ods, :material]
        setproperty!(par, item, InitParameters(item))
    end
    return par
end

function InitParameters(::Type{Val{:general}})
    general = InitParameters(nothing)
    general.casename = Entry(String, "", "Mnemonic name of the case being run")
    options = [
        :ods => "Load data from ODS saved in .json format",
        :scalars => "Initialize FUSE run form scalar FUSE parameters",
        :gasc => "Initialize FUSE run form GASC output file saved in .json format",
    ]
    general.init_from = Switch(options, "", "Initialize run from")
    return general
end

function InitParameters(::Type{Val{:material}})
    material = InitParameters(nothing)
    material.wall = Switch(FusionMaterials.available_materials("wall_materials"), "", "Material used for the wall"; default="Steel, Stainless 316")
    material.blanket = Switch(FusionMaterials.available_materials("blanket_materials"), "", "Material used for blanket coils")
    material.shield = Switch(FusionMaterials.available_materials("shield_materials"), "", "Material used for the shield")
    return material
end

function InitParameters(::Type{Val{:equilibrium}})
    equilibrium = InitParameters(nothing)
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

function InitParameters(::Type{Val{:core_profiles}})
    core_profiles = InitParameters(nothing)
    core_profiles.ne_ped = Entry(Real, "m^-3", "Pedestal electron density")
    core_profiles.n_peaking = Entry(Real, "", "Ratio of core/pedestal densities")
    core_profiles.T_shaping = Entry(Real, "", "Temperature shaping factor")
    core_profiles.w_ped = Entry(Real, "", "Pedestal width expressed in fraction of ψₙ")
    core_profiles.zeff = Entry(Real, "", "Effective ion charge")
    core_profiles.rot_core = Entry(Real, IMAS.core_profiles__profiles_1d, :rotation_frequency_tor_sonic)
    core_profiles.ngrid = Entry(Int, "", "Resolution of the core_profiles grid"; default=101)
    core_profiles.bulk = Entry(Symbol, "", "Bulk ion species")
    core_profiles.impurity = Entry(Symbol, "", "Impurity ion species")
    core_profiles.ejima = Entry(Real, "", "Ejima coefficient"; default=0.4)
    return core_profiles
end

function InitParameters(::Type{Val{:pf_active}})
    pf_active = InitParameters(nothing)
    pf_active.n_oh_coils = Entry(Int, "", "Number of OH coils")
    pf_active.n_pf_coils_inside = Entry(Int, "", "Number of PF coils inside of the TF")
    pf_active.n_pf_coils_outside = Entry(Int, "", "Number of PF coils outside of the TF")
    pf_active.technology = InitParameters(:coil_technology)
    return pf_active
end

function InitParameters(::Type{Val{:tf}})
    tf = InitParameters(nothing)
    tf.n_coils = Entry(Int, "", "Number of TF coils")
    options = [:princeton_D_exact, :princeton_D_approx, :princeton_D_scaled, :rectangle, :triple_arc, :miller, :spline]
    tf.shape = Switch(options, "", "Shape of the TF coils"; default=:princeton_D_scaled)
    tf.technology = InitParameters(:coil_technology)
    return tf
end

function InitParameters(::Type{Val{:oh}})
    oh = InitParameters(nothing)
    oh.technology = InitParameters(:coil_technology)
    oh.flattop_duration = Entry(Real, "s", "Duration of the flattop (use Inf for steady-state)")
    return oh
end

function InitParameters(::Type{Val{:center_stack}})
    center_stack = InitParameters(nothing)
    center_stack.bucked = Entry(Bool, "", "flag for bucked boundary conditions between TF and OH (and center plug, if present"; default=false)
    center_stack.noslip = Entry(Bool, "", "flag for no slip conditions between TF and OH (and center plug, if present)"; default=false)
    center_stack.plug = Entry(Bool, "", "flag for center plug"; default=false)
    return center_stack
end


function InitParameters(::Type{Val{:nbi}})
    nbi = InitParameters(nothing)
    nbi.power_launched = Entry(Union{X,Vector{X}} where {X<:Real}, "W", "Beam power")
    nbi.beam_energy = Entry(Union{X,Vector{X}} where {X<:Real}, "eV", "Beam energy")
    nbi.beam_mass = Entry(Union{X,Vector{X}} where {X<:Real}, "AU", "Beam mass"; default=2.0)
    nbi.toroidal_angle = Entry(Union{X,Vector{X}} where {X<:Real}, "rad", "toroidal angle of injection"; default=0.0)
    return nbi
end

function InitParameters(::Type{Val{:ec}})
    ec = InitParameters(nothing)
    ec.power_launched = Entry(Union{X,Vector{X}} where {X<:Real}, "W", "EC launched power")
    return ec
end

function InitParameters(::Type{Val{:ic}})
    ic = InitParameters(nothing)
    ic.power_launched = Entry(Union{X,Vector{X}} where {X<:Real}, "W", "IC launched power")
    return ic
end

function InitParameters(::Type{Val{:lh}})
    lh = InitParameters(nothing)
    lh.power_launched = Entry(Union{X,Vector{X}} where {X<:Real}, "W", "LH launched power")
    return lh
end

function InitParameters(::Type{Val{:build}})
    build = InitParameters(nothing)
    build.layers = Entry(DataStructures.OrderedDict, "m", "Sorted dictionary of layers thicknesses in radial build")
    build.blanket = Entry(Float64, "", "Fraction of blanket in radial build")
    build.shield = Entry(Float64, "", "Fraction of shield in radial build")
    build.vessel = Entry(Float64, "", "Fraction of vessel in radial build")
    build.symmetric = Entry(Bool, "", "Is the build up-down symmetric")
    return build
end

function InitParameters(::Type{Val{:gasc}})
    gasc = InitParameters(nothing)
    gasc.filename = Entry(String, "", "Output GASC .json file from which data will be loaded")
    gasc.case = Entry(Int, "", "Number of the GASC run to load")
    gasc.no_small_gaps = Entry(Bool, "", "Remove small gaps from the GASC radial build"; default=true)
    return gasc
end

function InitParameters(::Type{Val{:ods}})
    ods = InitParameters(nothing)
    ods.filename = Entry(String, "", "ODS.json file from which equilibrium is loaded")
    return ods
end

function InitParameters(::Type{Val{:coil_technology}})
    coil_tech = InitParameters(nothing)
    coil_tech.material = Switch(FusionMaterials.available_materials("magnet_materials"), "", "Technology used for the coil.")
    coil_tech.temperature = Entry(Real, "K", "Coil temperature")
    coil_tech.thermal_strain = Entry(Real, "", "Fraction of thermal expansion strain over maximum total strain on coil")
    coil_tech.JxB_strain = Entry(Real, "", "Fraction of maximum JxB strain over maximum total strain on coil")
    coil_tech.fraction_stainless = Entry(Real, "", "Fraction of stainless steel in the coil cross-sectional areas")
    coil_tech.fraction_void = Entry(Real, "", "Fraction of `void` in the coil cross-sectional area. Void is everything (like coolant) that is not structural nor conductor.")
    coil_tech.ratio_SC_to_copper = Entry(Real, "", "Fraction of superconductor to copper cross-sectional areas")
    return coil_tech
end

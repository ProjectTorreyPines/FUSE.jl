using FusionMaterials: FusionMaterials
import SimulationParameters: SwitchOption

import IMAS: BuildLayerType, _plasma_, _gap_, _oh_, _tf_, _shield_, _blanket_, _wall_, _vessel_, _cryostat_, _divertor_, _port_
import IMAS: BuildLayerSide, _lfs_, _lhfs_, _hfs_, _in_, _out_
import IMAS: BuildLayerShape, _offset_, _negative_offset_, _convex_hull_, _princeton_D_exact_, _princeton_D_, _princeton_D_scaled_, _rectangle_, _double_ellipse_, _circle_ellipse_,
    _triple_arc_, _miller_, _silo_, _racetrack_, _undefined_

const layer_shape_options = Dict(Symbol(string(e)[2:end-1]) => SwitchOption(e, string(e)[2:end-1]) for e in instances(IMAS.BuildLayerShape))
const layer_type_options = Dict(Symbol(string(e)[2:end-1]) => SwitchOption(e, string(e)[2:end-1]) for e in instances(IMAS.BuildLayerType))
const layer_side_options = Dict(Symbol(string(e)[2:end-1]) => SwitchOption(e, string(e)[2:end-1]) for e in instances(IMAS.BuildLayerSide))

Base.@kwdef mutable struct FUSEparameters__time{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :time
    pulse_shedule_time_basis::Entry{AbstractRange{Float64}} = Entry{AbstractRange{Float64}}("s", "Time basis used to discretize the pulse schedule")
    simulation_start::Entry{Float64} = Entry{Float64}("s", "Time at which the simulation starts"; default=0.0)
end

Base.@kwdef mutable struct FUSEparameters__general{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :general
    casename::Entry{String} = Entry{String}("-", "Sort mnemonic name of the case being run")
    description::Entry{String} = Entry{String}("-", "Longer description of the case being run")
    init_from::Switch{Symbol} = Switch{Symbol}(
        [
            :ods => "Load data from ODS saved in .json format (where possible, and fallback on scalars otherwise)",
            :scalars => "Initialize FUSE run from scalar parameters"
        ], "-", "Initialize run from")
    dd::Entry{IMAS.dd} = Entry{IMAS.dd}("-", "`dd` to initialize from")
end

Base.@kwdef mutable struct FUSEparameters__rampup{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :rampup
    side::Switch{Symbol} = Switch{Symbol}([:hfs, :lfs], "-", "Side of the vacuum vessel where the plasma is limited at breakdown")
    ends_at::Entry{Float64} = Entry{Float64}("s", "Until when does the rampup lasts")
    diverted_at::Entry{Float64} = Entry{Float64}("s", "Time at which x-point is formed and plasma can peel-off the wall")
end

Base.@kwdef mutable struct FUSEparameters__equilibrium{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :equilibrium
    B0::Entry{T} = Entry{T}(IMAS.equilibrium__vacuum_toroidal_field, :b0)
    R0::Entry{T} = Entry{T}("m", "Geometric genter of the plasma. NOTE: This also scales the radial build layers."; check=x -> @assert x > 0.0 "must be: R0 > 0.0")
    Z0::Entry{T} = Entry{T}("m", "Z offset of the machine midplane"; default=0.0)
    ϵ::Entry{T} = Entry{T}("-", "Plasma inverse aspect ratio (a/R0). NOTE: This also scales the radial build layers."; check=x -> @assert 0.0 < x < 1.0 "must be: 0.0 < ϵ < 1.0")
    κ::Entry{T} = Entry{T}("-", "Plasma elongation. NOTE: If < 1.0 it defines the fraction of maximum controllable elongation estimate.")
    tilt::Entry{T} = Entry{T}("-", "Tilt of the plasma boundary [MXH c0]"; default=0.0)
    δ::Entry{T} = Entry{T}("-", "Triangularity of the plasma boundary [MXH sin(s1)]"; default=0.0)
    ζ::Entry{T} = Entry{T}("-", "Squareness of the plasma boundary [MXH -s2]"; default=0.0)
    𝚶::Entry{T} = Entry{T}("-", "Ovality of the plasma boundary [MXH c1]"; default=0.0)
    twist::Entry{T} = Entry{T}("-", "Twist of the plasma boundary [MXH c2]"; default=0.0)
    pressure_core::Entry{T} = Entry{T}("Pa", "On axis pressure"; check=x -> @assert x > 0.0 "must be: P > 0.0")
    ip::Entry{T} = Entry{T}(IMAS.equilibrium__time_slice___global_quantities, :ip)
    xpoints::Switch{Symbol} = Switch{Symbol}([:lower, :upper, :double, :none], "-", "X-points configuration")
    ngrid::Entry{Int} = Entry{Int}("-", "Resolution of the equilibrium grid"; default=129)
    field_null_surface::Entry{T} =
        Entry{T}("-", "ψn value of the field_null_surface. Disable with 0.0"; default=0.75, check=x -> @assert x > 0.0 "must be: field_null_surface > 0.0")
    boundary_from::Switch{Symbol} = Switch{Symbol}([:scalars, :MXH_params, :rz_points, :ods], "-", "The starting r, z boundary taken from")
    MXH_params::Entry{Vector{T}} = Entry{Vector{T}}("-", "Vector of MXH flats")
    rz_points::Entry{Vector{Vector{T}}} = Entry{Vector{Vector{T}}}("m", "R_Z boundary as Vector{Vector{$T}}} : r = rz_points[1], z = rz_points[2]")
end

Base.@kwdef mutable struct FUSEparameters__core_profiles{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :core_profiles
    plasma_mode::Switch{Symbol} = Switch{Symbol}([:H_mode, :L_mode], "-", "Plasma configuration"; default=:H_mode)
    ne_value::Entry{T} = Entry{T}("-", "Value based on setup method"; check=x -> @assert x > 0.0 "must be > 0.0")
    ne_setting::Switch{Symbol} = Switch{Symbol}([:ne_ped, :ne_line, :greenwald_fraction, :greenwald_fraction_ped], "-", "Way to set the electron density")
    w_ped::Entry{T} = Entry{T}("-", "Pedestal width expressed in fraction of ψₙ"; default=0.05, check=x -> @assert x > 0.0 "must be: w_ped > 0.0")
    ne_sep_to_ped_ratio::Entry{T} =
        Entry{T}("-", "Ratio used to set the sepeartrix density based on the pedestal density"; default=0.25, check=x -> @assert x > 0.0 "must be: ne_sep_to_ped_ratio > 0.0")
    T_ratio::Entry{T} = Entry{T}("-", "Ti/Te ratio"; check=x -> @assert x > 0.0 "must be: T_ratio > 0.0")
    T_shaping::Entry{T} = Entry{T}("-", "Temperature shaping factor")
    n_shaping::Entry{T} = Entry{T}("-", "Density shaping factor")
    zeff::Entry{T} = Entry{T}("-", "Effective ion charge"; check=x -> @assert x >= 1.0 "must be: zeff > 1.0")
    rot_core::Entry{T} = Entry{T}(IMAS.core_profiles__profiles_1d, :rotation_frequency_tor_sonic)
    ngrid::Entry{Int} = Entry{Int}("-", "Resolution of the core_profiles grid"; default=101, check=x -> @assert x >= 11 "must be: ngrid >= 11")
    bulk::Entry{Symbol} = Entry{Symbol}("-", "Bulk ion species")
    impurity::Entry{Symbol} = Entry{Symbol}("-", "Impurity ion species")
    helium_fraction::Entry{T} = Entry{T}("-", "Helium density / electron density fraction"; check=x -> @assert 0.0 <= x <= 0.5 "must be: 0.0 <= helium_fraction <= 0.5")
    ejima::Entry{T} = Entry{T}("-", "Ejima coefficient"; default=0.4, check=x -> @assert 0.0 <= x < 1.0 "must be: 0.0 <= ejima < 1.0")
    polarized_fuel_fraction::Entry{T} = Entry{T}("-", "Spin polarized fuel fraction"; default=0.0, check=x -> @assert 0.0 < x < 1.0 "must be: 0.0 < polarized_fuel_fraction < 1.0")
end

Base.@kwdef mutable struct FUSEparameters__pf_active{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :pf_active
    n_coils_inside::Entry{Int} = Entry{Int}("-", "Number of PF coils inside of the TF"; check=x -> @assert x >= 0 "must be: n_coils_inside >= 0")
    n_coils_outside::Entry{Int} = Entry{Int}("-", "Number of PF coils outside of the TF"; check=x -> @assert x >= 0 "must be: n_coils_outside >= 0")
    technology::Switch{Symbol} = Switch{Symbol}(FusionMaterials.supported_coil_techs(), "-", "PF coils technology")
end

Base.@kwdef mutable struct FUSEparameters__nb_unit{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :nbi
    power_launched::Entry{T} = Entry{T}("W", "Beam power"; check=x -> @assert x >= 0.0 "must be: power_launched >= 0.0")
    rho_0::Entry{T} = Entry{T}("-", "Desired radial location of the deposition profile"; default=0.0, check=x -> @assert x >= 0.0 "must be: rho_0 >= 0.0")
    width::Entry{T} = Entry{T}("-", "Desired width of the deposition profile"; default=0.3, check=x -> @assert x > 0.0 "must be: width > 0.0")
    beam_energy::Entry{T} = Entry{T}("eV", "Beam energy"; check=x -> @assert x >= 0.0 "must be: beam_energy >= 0.0")
    beam_mass::Entry{T} = Entry{T}("AU", "Beam mass"; default=2.0, check=x -> @assert x >= 1.0 "must be: beam_mass >= 1.0")
    toroidal_angle::Entry{T} = Entry{T}("rad", "Toroidal angle of injection"; default=0.0)
    efficiency_conversion::Entry{T} = Entry{T}(IMAS.nbi__unit___efficiency, :conversion; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_conversion > 0.0")
    efficiency_transmission::Entry{T} = Entry{T}(IMAS.nbi__unit___efficiency, :transmission; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_transmission > 0.0")
end

Base.@kwdef mutable struct FUSEparameters__ec_launcher{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :ec_launcher
    power_launched::Entry{T} = Entry{T}("W", "EC launched power"; check=x -> @assert x >= 0.0 "must be: power_launched >= 0.0")
    rho_0::Entry{T} = Entry{T}("-", "Desired radial location of the deposition profile"; default=0.5, check=x -> @assert x >= 0.0 "must be: rho_0 >= 0.0")
    width::Entry{T} = Entry{T}("-", "Desired width of the deposition profile"; default=0.025, check=x -> @assert x > 0.0 "must be: width > 0.0")
    efficiency_conversion::Entry{T} = Entry{T}(IMAS.ec_launchers__beam___efficiency, :conversion; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_conversion > 0.0")
    efficiency_transmission::Entry{T} =
        Entry{T}(IMAS.ec_launchers__beam___efficiency, :transmission; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_transmission > 0.0")
end

Base.@kwdef mutable struct FUSEparameters__ic_antenna{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :ic_antenna
    power_launched::Entry{T} = Entry{T}("W", "IC launched power"; check=x -> @assert x >= 0.0 "must be: power_launched >= 0.0")
    rho_0::Entry{T} = Entry{T}("-", "Desired radial location of the deposition profile"; default=0.0, check=x -> @assert x >= 0.0 "must be: rho_0 >= 0.0")
    width::Entry{T} = Entry{T}("-", "Desired width of the deposition profile"; default=0.1, check=x -> @assert x > 0.0 "must be: width > 0.0")
    efficiency_conversion::Entry{T} = Entry{T}(IMAS.ic_antennas__antenna___efficiency, :conversion; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_conversion > 0.0")
    efficiency_transmission::Entry{T} =
        Entry{T}(IMAS.ic_antennas__antenna___efficiency, :transmission; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_transmission > 0.0")
    efficiency_coupling::Entry{T} = Entry{T}(IMAS.ic_antennas__antenna___efficiency, :coupling; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_coupling > 0.0")
end

Base.@kwdef mutable struct FUSEparameters__lh_antenna{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :lh_antenna
    power_launched::Entry{T} = Entry{T}("W", "LH launched power"; check=x -> @assert x >= 0.0 "must be: power_launched >= 0.0")
    rho_0::Entry{T} = Entry{T}("-", "Desired radial location of the deposition profile"; default=0.8, check=x -> @assert x >= 0.0 "must be: rho_0 >= 0.0")
    width::Entry{T} = Entry{T}("-", "Desired width of the deposition profile"; default=0.05, check=x -> @assert x > 0.0 "must be: width > 0.0")
    efficiency_conversion::Entry{T} = Entry{T}(IMAS.lh_antennas__antenna___efficiency, :conversion; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_conversion > 0.0")
    efficiency_transmission::Entry{T} =
        Entry{T}(IMAS.lh_antennas__antenna___efficiency, :transmission; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_transmission > 0.0")
    efficiency_coupling::Entry{T} = Entry{T}(IMAS.lh_antennas__antenna___efficiency, :coupling; default=1.0)
end

Base.@kwdef mutable struct FUSEparameters__pellet_launcher{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :pellet_launcher
    frequency::Entry{T} = Entry{T}("Hz", "Frequency of pellets launched"; check=x -> @assert x >= 0.0 "pellet frequency must be >= 0.0")
    rho_0::Entry{T} = Entry{T}("-", "Desired radial location of the deposition profile"; default=0.5, check=x -> @assert x >= 0.0 "must be: rho_0 >= 0.0")
    width::Entry{T} = Entry{T}("-", "Desired width of the deposition profile"; default=0.25, check=x -> @assert x > 0.0 "must be: width > 0.0")
    shape::Switch{Symbol} = Switch{Symbol}([:spherical, :cylindrical, :rectangular], "-", "The pellet geometry"; default=:spherical)
    species::Switch{Symbol} = Switch{Symbol}([:H, :D, :T, :DT, :C, :Ne], "-", "Pellet species")
    size::Entry{Vector{T}} = Entry{Vector{T}}(
        "m",
        "Vector of geometric dimensions describing the pellet size for a given shape (spherical: [r], cylindrical: [d, l], rectangular: [x,y,z])";
        check=x -> @assert all(x .> 0.0) "All pellet shape dimensions must be > 0.0"
    )
end

Base.@kwdef mutable struct FUSEparameters__ods{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :ods
    filename::Entry{String} = Entry{String}("-", "ODS.json file(s) from which equilibrium is loaded. Multiple comma-separated ODSs can be specified.")
end

Base.@kwdef mutable struct FUSEparameters__tf{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :tf
    n_coils::Entry{Int} = Entry{Int}("-", "Number of TF coils"; check=x -> @assert x >= 0 "must be: n_coils >= 0")
    shape::Switch{BuildLayerShape} = Switch{BuildLayerShape}(layer_shape_options, "-", "Shape of the TF coils")
    ripple::Entry{T} =
        Entry{T}("-", "Fraction of toroidal field ripple evaluated at the outermost radius of the plasma chamber"; default=0.01, check=x -> @assert x > 0.0 "must be: ripple > 0.0")
    technology::Switch{Symbol} = Switch{Symbol}(FusionMaterials.supported_coil_techs(), "-", "TF coils technology")
end

Base.@kwdef mutable struct FUSEparameters__oh{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :oh
    n_coils::Entry{Int} = Entry{Int}("-", "Number of OH coils"; check=x -> @assert x >= 0 "must be: n_coils >= 0")
    technology::Switch{Symbol} = Switch{Symbol}(FusionMaterials.supported_coil_techs(), "-", "OH coils technology")
end

Base.@kwdef mutable struct FUSEparameters__center_stack{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :center_stack
    bucked::Entry{Bool} = Entry{Bool}("-", "Flag for bucked boundary conditions between TF and OH (and center plug, if present)"; default=false)
    noslip::Entry{Bool} = Entry{Bool}("-", "Flag for no slip conditions between TF and OH (and center plug, if present)"; default=false)
    plug::Entry{Bool} = Entry{Bool}("-", "Flag for center plug"; default=false)
end

Base.@kwdef mutable struct FUSEparameters__build_layer{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :layer
    name::Entry{String} = Entry{String}("-", "Name of the layer")
    thickness::Entry{Float64} =
        Entry{Float64}("-", "Relative thickness of the layer (layers actual thickness is scaled to match plasma R0)"; check=x -> @assert x >= 0.0 "must be: thickness >= 0.0")
    material::Switch{Symbol} = Switch{Symbol}(
        FusionMaterials.all_materials(),
        "-",
        "Material of the layer"
    )
    shape::Switch{BuildLayerShape} = Switch{BuildLayerShape}(layer_shape_options, "-", "Shape of the layer")
    type::Switch{BuildLayerType} = Switch{BuildLayerType}(layer_type_options, "-", "Type of the layer")
    side::Switch{BuildLayerSide} = Switch{BuildLayerSide}(layer_side_options, "-", "Side of the layer")
end

Base.@kwdef mutable struct FUSEparameters__build{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :build
    layers::ParametersVector{FUSEparameters__build_layer{T}} = ParametersVector{FUSEparameters__build_layer{T}}()
    plasma_gap::Entry{T} =
        Entry{T}("-", "Fraction of vacuum gap between first wall and plasma separatrix in radial build"; default=0.1, check=x -> @assert x > 0.0 "must be: plasma_gap > 0.0")
    symmetric::Entry{Bool} = Entry{Bool}("-", "Is the build up-down symmetric")
    divertors::Switch{Symbol} = Switch{Symbol}([:lower, :upper, :double, :none, :from_x_points], "-", "Divertors configuration"; default=:from_x_points)
    n_first_wall_conformal_layers::Entry{Int} =
        Entry{Int}("-", "Number of layers that are conformal to the first wall"; default=1, check=x -> @assert x > 0 "must be: n_first_wall_conformal_layers > 0")
end

Base.@kwdef mutable struct FUSEparameters__requirements{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :requirements
    power_electric_net::Entry{T} = Entry{T}(IMAS.requirements, :power_electric_net; check=x -> @assert x >= 0.0 "must be: power_electric_net >= 0.0")
    flattop_duration::Entry{T} = Entry{T}(IMAS.requirements, :flattop_duration; check=x -> @assert x >= 0.0 "must be: flattop_duration >= 0.0")
    log10_flattop_duration::Entry{T} =
        Entry{T}("log10(s)", "Log10 value of the duration of the flattop (use Inf for steady-state). Preferred over `flattop_duration` for optimization studies.")
    tritium_breeding_ratio::Entry{T} = Entry{T}(IMAS.requirements, :tritium_breeding_ratio; check=x -> @assert x >= 0.0 "must be: tritium_breeding_ratio >= 0.0")
    cost::Entry{T} = Entry{T}(IMAS.requirements, :cost; check=x -> @assert x >= 0.0 "must be: cost >= 0.0")
    ne_peaking::Entry{T} = Entry{T}(IMAS.requirements, :ne_peaking; check=x -> @assert x >= 0.0 "must be: ne_peaking >= 0.0")
    q_pol_omp::Entry{T} = Entry{T}(IMAS.requirements, :q_pol_omp; check=x -> @assert x >= 0.0 "must be: q_pol_omp >= 0.0")
    lh_power_threshold_fraction::Entry{T} = Entry{T}(IMAS.requirements, :lh_power_threshold_fraction; check=x -> @assert x >= 0.0 "must be: lh_power_threshold_fraction >= 0.0")
    h98y2::Entry{T} = Entry{T}(IMAS.requirements, :h98y2; check=x -> @assert x >= 0.0 "must be: h98y2 >= 0.0")
    hds03::Entry{T} = Entry{T}(IMAS.requirements, :hds03; check=x -> @assert x >= 0.0 "must be: hds03 >= 0.0")
    beta_normal::Entry{T} = Entry{T}(IMAS.requirements, :beta_normal; check=x -> @assert x >= 0.0 "must be: βn >= 0.0")
    Psol_R::Entry{T} = Entry{T}(IMAS.requirements, :Psol_R; check=x -> @assert x >= 0.0 "must be: Psol/R >= 0.0")
    q95::Entry{T} = Entry{T}(IMAS.requirements, :q95; check=x -> @assert x >= 0.0 "must be: q95 >= 0.0")
    coil_j_margin::Entry{T} = Entry{T}(IMAS.requirements, :coil_j_margin; default=0.4, check=x -> @assert x >= 0.0 "must be: coil_j_margin >= 0.0")
    coil_stress_margin::Entry{T} = Entry{T}(IMAS.requirements, :coil_stress_margin; default=0.2, check=x -> @assert x >= 0.0 "must be: coil_j_margin >= 0.0")
end

Base.@kwdef mutable struct FUSEparameters__balance_of_plant{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :balance_of_plant
    cycle_type::Switch{Symbol} = Switch{Symbol}([:rankine, :brayton], "-", "Thermal cycle type"; default=:rankine)
end

mutable struct ParametersInits{T<:Real} <: ParametersAllInits{T}
    _parent::WeakRef
    _name::Symbol
    general::FUSEparameters__general{T}
    time::FUSEparameters__time{T}
    ods::FUSEparameters__ods{T}
    equilibrium::FUSEparameters__equilibrium{T}
    core_profiles::FUSEparameters__core_profiles{T}
    pf_active::FUSEparameters__pf_active{T}
    rampup::FUSEparameters__rampup{T}
    nb_unit::ParametersVector{FUSEparameters__nb_unit{T}}
    ec_launcher::ParametersVector{FUSEparameters__ec_launcher{T}}
    pellet_launcher::ParametersVector{FUSEparameters__pellet_launcher{T}}
    ic_antenna::ParametersVector{FUSEparameters__ic_antenna{T}}
    lh_antenna::ParametersVector{FUSEparameters__lh_antenna{T}}
    build::FUSEparameters__build{T}
    center_stack::FUSEparameters__center_stack{T} #
    tf::FUSEparameters__tf{T}
    oh::FUSEparameters__oh{T}
    bop::FUSEparameters__balance_of_plant{T}
    requirements::FUSEparameters__requirements{T}
end

function ParametersInits{T}(; n_nb::Int=0, n_ec::Int=0, n_pl::Int=0, n_ic::Int=0, n_lh::Int=0, n_layers::Int=0) where {T<:Real}
    ini = ParametersInits{T}(
        WeakRef(nothing),
        :ini,
        FUSEparameters__general{T}(),
        FUSEparameters__time{T}(),
        FUSEparameters__ods{T}(),
        FUSEparameters__equilibrium{T}(),
        FUSEparameters__core_profiles{T}(),
        FUSEparameters__pf_active{T}(),
        FUSEparameters__rampup{T}(),
        ParametersVector{FUSEparameters__nb_unit{T}}(),
        ParametersVector{FUSEparameters__ec_launcher{T}}(),
        ParametersVector{FUSEparameters__pellet_launcher{T}}(),
        ParametersVector{FUSEparameters__ic_antenna{T}}(),
        ParametersVector{FUSEparameters__lh_antenna{T}}(),
        FUSEparameters__build{T}(),
        FUSEparameters__center_stack{T}(),
        FUSEparameters__tf{T}(),
        FUSEparameters__oh{T}(),
        FUSEparameters__balance_of_plant{T}(),
        FUSEparameters__requirements{T}())

    for k in 1:n_layers
        push!(ini.build.layers, FUSEparameters__build_layer{T}())
    end

    for k in 1:n_nb
        push!(ini.nb_unit, FUSEparameters__nb_unit{T}())
    end

    for k in 1:n_ec
        push!(ini.ec_launcher, FUSEparameters__ec_launcher{T}())
    end

    for k in 1:n_pl
        push!(ini.pellet_launcher, FUSEparameters__pellet_launcher{T}())
    end

    for k in 1:n_ic
        push!(ini.ic_antenna, FUSEparameters__ic_antenna{T}())
    end

    for k in 1:n_lh
        push!(ini.lh_antenna, FUSEparameters__lh_antenna{T}())
    end

    setup_parameters!(ini)

    return ini
end

function ParametersInits(; kw...)
    return ParametersInits{Float64}(; kw...)
end

################################
# functions for populating ini #
################################
"""
    ini_equilibrium_elongation_true(equilibrium::FUSEparameters__equilibrium)

if elongation <1.0 then expresses elongation as fraction of maximum controllable elongation estimate
"""
function ini_equilibrium_elongation_true(equilibrium::FUSEparameters__equilibrium)
    if !ismissing(equilibrium, :κ)
        if equilibrium.κ < 1.0 && !ismissing(equilibrium, :ϵ)
            return IMAS.elongation_limit(1.0 / equilibrium.ϵ) * equilibrium.κ
        else
            return equilibrium.κ
        end
    else
        return missing
    end
end

"""
    ini_equilibrium_elongation_true(κ::T, ϵ::T) where {T<:Real}

if elongation <1.0 then expresses elongation as fraction of maximum controllable elongation estimate
"""
function ini_equilibrium_elongation_true(κ::T, ϵ::T) where {T<:Real}
    if κ < 1.0
        return IMAS.elongation_limit(1.0 / ϵ) * κ
    else
        return κ
    end
end

"""
    (equilibrium::FUSEparameters__equilibrium)(mxh::IMAS.MXH)

ini.equilibrium scalars from MXH parametrization
"""
function (equilibrium::FUSEparameters__equilibrium)(mxh::IMAS.MXH)
    equilibrium.ϵ = mxh.ϵ
    equilibrium.R0 = mxh.R0
    equilibrium.Z0 = mxh.Z0
    equilibrium.κ = mxh.κ
    equilibrium.tilt = mxh.c0
    if length(mxh.s) >= 1
        equilibrium.δ = sin(mxh.s[1])
    end
    if length(mxh.s) >= 2
        equilibrium.ζ = -mxh.s[2]
    end
    if length(mxh.c) >= 1
        equilibrium.𝚶 = mxh.c[1]
    end
    if length(mxh.c) >= 2
        equilibrium.twist = mxh.c[2]
    end
    return equilibrium
end

"""
    MXHboundary(ini::ParametersAllInits)::MXHboundary

return MHXboundary representation of independenty of how it was input in ini.equilibrium
"""
function MXHboundary(ini::ParametersAllInits; kw...)::MXHboundary
    if ini.general.init_from == :ods
        dd = load_ods(ini)
    else
        dd = IMAS.dd()
    end
    return MXHboundary(ini, dd; kw...)
end

function MXHboundary(ini::ParametersAllInits, dd::IMAS.dd; kw...)::MXHboundary
    init_from = ini.general.init_from
    if init_from == :ods
        if !ismissing(dd.equilibrium, :time) && length(dd.equilibrium.time) > 0
            dd.global_time = ini.time.simulation_start
            eqt = dd.equilibrium.time_slice[]
            fw = IMAS.first_wall(dd.wall)
            IMAS.flux_surfaces(eqt, fw.r, fw.z)
        else
            init_from = :scalars
        end
    end

    boundary_from = ini.equilibrium.boundary_from
    if boundary_from == :ods
        pr, pz = eqt.boundary.outline.r, eqt.boundary.outline.z
        pr, pz = IMAS.resample_plasma_boundary(pr, pz; n_points=101)
        pr, pz = IMAS.reorder_flux_surface!(pr, pz)
        mxh = IMAS.MXH(pr, pz, 4)

    elseif boundary_from == :rz_points
        # R,Z boundary from points
        if ismissing(ini.equilibrium, :rz_points)
            error("ini.equilibrium.boundary_from is set as $boundary_from but rz_points wasn't set")
        end
        pr, pz = ini.equilibrium.rz_points[1], ini.equilibrium.rz_points[2]
        pr, pz = IMAS.resample_plasma_boundary(pr, pz; n_points=101)
        pr, pz = IMAS.reorder_flux_surface!(pr, pz)
        mxh = IMAS.MXH(pr, pz, 4)

    elseif boundary_from == :MXH_params
        # R,Z boundary from MXH
        if ismissing(ini.equilibrium, :MXH_params)
            error("ini.equilibrium.boundary_from is set as $boundary_from but MXH_params wasn't set")
        end
        mxh = IMAS.MXH(ini.equilibrium.MXH_params)

    elseif boundary_from == :scalars
        # R,Z boundary from scalars
        mxh = IMAS.MXH(
            ini.equilibrium.R0,
            ini.equilibrium.Z0,
            ini.equilibrium.ϵ,
            ini_equilibrium_elongation_true(ini.equilibrium),
            ini.equilibrium.tilt,
            [ini.equilibrium.𝚶, 0.0],
            [asin(ini.equilibrium.δ), -ini.equilibrium.ζ])
    else
        error("ini.equilibrium.boundary_from must be one of [:scalars, :rz_points, :MXH_params, :ods]")
    end

    if boundary_from == :ods
        # in case of ODS we have all information to generate MXHboundary
        RX = [x_point.r for x_point in eqt.boundary.x_point]
        ZX = [x_point.z for x_point in eqt.boundary.x_point]
        mxhb = MXHboundary(mxh, RX, ZX, pr, pz)
    else
        # all other cases we must reconcile mxh boundary with requested x-points
        nx = n_xpoints(ini.equilibrium.xpoints)
        mxhb = fitMXHboundary(mxh, nx; kw...)
    end

    return mxhb
end

function n_xpoints(xpoints::Symbol)
    if xpoints == :none
        return 0
    elseif xpoints == :lower
        return -1
    elseif xpoints == :upper
        return 1
    elseif xpoints == :double
        return 2
    else
        error("xpoints can only be [:none, :lower, :upper, :double]")
    end
end

"""
    load_ods(ini::ParametersAllInits)

Load ODSs as specified in `ini.ods.filename`
and sets `dd.global_time` equal to `ini.time.simulation_start`

NOTE: supports multiple comma-separated filenames
"""
function load_ods(ini::ParametersAllInits)
    dd = load_ods(ini.ods.filename)
    dd.global_time = ini.time.simulation_start
    for field in keys(dd)
        ids = getproperty(dd, field)
        if !ismissing(ids, :time) && length(ids.time) == 1
            IMAS.retime!(ids, dd.global_time)
        end
    end
    return dd
end

"""
    load_ods(filenames::String)

Load multiple comma-separated filenames into a single dd
"""
function load_ods(filenames::String)
    return load_ods(strip.(split(filenames, ",")))
end

"""
    load_ods(filenames::Vector{String})

Load multiple ODSs into a single `dd`
"""
function load_ods(filenames::Vector{<:AbstractString})
    dd = IMAS.dd()
    for filename in filenames
        filename = replace(filename, r"^__FUSE__" => __FUSE__)
        dd1 = IMAS.json2imas(filename)
        merge!(dd, dd1)
    end
    IMAS.last_global_time(dd)
    return dd
end

###############
# save / load #
###############
"""
    ini2json(ini::ParametersAllInits, filename::AbstractString; kw...)

Save the FUSE parameters to a JSON file with give `filename`

`kw` arguments are passed to the JSON.print function
"""
function ini2json(ini::ParametersAllInits, filename::AbstractString; kw...)
    return SimulationParameters.par2json(ini, filename; kw...)
end

"""
    json2ini(filename::AbstractString, ini::ParametersAllInits=ParametersInits())

Load the FUSE parameters from a JSON file with given `filename`
"""
function json2ini(filename::AbstractString, ini::ParametersAllInits=ParametersInits())
    return SimulationParameters.json2par(filename, ini)
end

"""
    ini2yaml(ini::ParametersAllInits, filename::AbstractString; kw...)

Save the FUSE parameters to a YAML file with given `filename`

`kw` arguments are passed to the YAML.print function
"""
function ini2yaml(ini::ParametersAllInits, filename::AbstractString; kw...)
    return SimulationParameters.par2yaml(ini, filename; kw...)
end

"""
    yaml2ini(filename::AbstractString, ini::ParametersAllInits=ParametersInits())

Load the FUSE parameters from a YAML file with given `filename`
"""
function yaml2ini(filename::AbstractString, ini::ParametersAllInits=ParametersInits())
    return SimulationParameters.yaml2par(filename, ini)
end

########
# plot #
########
"""
    plot_ini(ini::ParametersAllInits; time0=global_time(ini))

Plots ini time dependent time traces including plasma boundary
"""
@recipe function plot_ini(ini::ParametersAllInits; time0=global_time(ini))
    @assert typeof(time0) <: Float64

    # count number of time-dependent parameters
    N = 0
    for par in SimulationParameters.leaves(ini)
        if typeof(par.value) <: Function
            N += 1
        end
    end
    layout := @layout [N + 1]

    time_bkp = ini.time.simulation_start
    try
        ini.time.simulation_start = time0

        # plot equilibrium including x-points
        mxhb = MXHboundary(ini)
        wr = wall_radii(mxhb.mxh.R0, mxhb.mxh.ϵ * mxhb.mxh.R0, ini.build.plasma_gap)
        @series begin
            label := ""
            seriestype := :vline
            subplot := 1
            colour := :black
            [wr.r_hfs, wr.r_lfs]
        end
        @series begin
            label := ""
            subplot := 1
            aspectratio := :equal
            xlim := (wr[1] - (wr[2] - wr[1]) / 10, wr[2] + (wr[2] - wr[1]) / 10)
            mxhb
        end

        # plot time dependent parameters
        k = 1
        for par in SimulationParameters.leaves(ini)
            if typeof(par.value) <: Function
                k += 1
                @series begin
                    label := ""
                    subplot := k
                    time0 := time0
                    par
                end
            end
        end

    finally
        ini.time.simulation_start = time_bkp
    end
end

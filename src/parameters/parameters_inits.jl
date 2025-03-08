using FusionMaterials: FusionMaterials
import SimulationParameters: SwitchOption

import IMAS: BuildLayerType, _plasma_, _gap_, _oh_, _tf_, _shield_, _blanket_, _wall_, _vessel_, _cryostat_, _divertor_, _port_
import IMAS: BuildLayerSide, _lfs_, _lhfs_, _hfs_, _in_, _out_
import IMAS: BuildLayerShape, _offset_, _negative_offset_, _convex_hull_, _princeton_D_, _mirror_princeton_D_, _princeton_D_scaled_, _mirror_princeton_D_scaled_, _rectangle_,
    _double_ellipse_, _mirror_double_ellipse_, _rectangle_ellipse_, _mirror_rectangle_ellipse_, _circle_ellipse_, _mirror_circle_ellipse_, _triple_arc_, _mirror_triple_arc_,
    _miller_, _silo_, _racetrack_, _undefined_

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
    pressure_core::Entry{T} =
        Entry{T}("Pa", "On axis pressure (NOTE: `pressure_core` can be calculated from `ini.core_profiles.Te_core`)"; check=x -> @assert x > 0.0 "must be: P > 0.0")
    ip::Entry{T} = Entry{T}(IMAS.equilibrium__time_slice___global_quantities, :ip)
    xpoints::Switch{Symbol} = Switch{Symbol}([:lower, :upper, :double, :none], "-", "X-points configuration")
    ngrid::Entry{Int} = Entry{Int}("-", "Resolution of the equilibrium grid"; default=129)
    field_null_surface::Entry{T} =
        Entry{T}("-", "ψn value of the field_null_surface. Disable with 0.0"; default=0.0, check=x -> @assert x >= 0.0 "must be: field_null_surface >= 0.0")
    boundary_from::Switch{Symbol} = Switch{Symbol}([:scalars, :MXH_params, :rz_points, :ods], "-", "The starting r, z boundary taken from")
    MXH_params::Entry{Vector{T}} = Entry{Vector{T}}("-", "Vector of MXH flats")
    rz_points::Entry{Vector{Vector{T}}} = Entry{Vector{Vector{T}}}("m", "R_Z boundary as Vector{Vector{$T}}} : r = rz_points[1], z = rz_points[2]")
end

Base.@kwdef mutable struct FUSEparameters__ITB{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :pf_active
    radius::Entry{T} = Entry{T}("-", "Rho at which the ITB starts")
    ne_width::Entry{T} = Entry{T}("-", "Width of the electron density ITB")
    ne_height_ratio::Entry{T} = Entry{T}("-", "Height of the electron density ITB, expressed as the ratio of the density without ITB evaluated on axis")
    Te_width::Entry{T} = Entry{T}("-", "Width of the electron temperature ITB")
    Te_height_ratio::Entry{T} = Entry{T}("-", "Height of the electron temperature ITB, expressed as the ratio of the temperature without ITB evaluated on axis")
end

Base.@kwdef mutable struct FUSEparameters__core_profiles{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :core_profiles
    plasma_mode::Switch{Symbol} = Switch{Symbol}([:H_mode, :L_mode], "-", "Plasma configuration"; default=:H_mode)
    w_ped::Entry{T} = Entry{T}("-", "Pedestal width expressed in fraction of ψₙ"; default=0.05, check=x -> @assert x > 0.0 "must be: w_ped > 0.0")
    ne_value::Entry{T} = Entry{T}("-", "Value based on setup method"; check=x -> @assert x > 0.0 "must be > 0.0")
    ne_setting::Switch{Symbol} = Switch{Symbol}([:ne_ped, :ne_line, :greenwald_fraction, :greenwald_fraction_ped], "-", "Way to set the electron density")
    ne_sep_to_ped_ratio::Entry{T} =
        Entry{T}(
            "-",
            "Ratio used to set the sepeartrix density based on the pedestal density";
            default=0.25,
            check=x -> @assert 1.0 > x > 0.0 "must be: 1.0 > ne_sep_to_ped_ratio > 0.0"
        )
    ne_core_to_ped_ratio::Entry{T} =
        Entry{T}("-", "Ratio used to set the core density based on the pedestal density"; default=1.4, check=x -> @assert x > 0.0 "must be: ne_core_to_ped_ratio > 0.0")
    ne_shaping::Entry{T} = Entry{T}("-", "Density shaping factor")
    Ti_Te_ratio::Entry{T} = Entry{T}("-", "Ti/Te ratio"; check=x -> @assert x > 0.0 "must be: Ti_Te_ratio > 0.0")
    Te_shaping::Entry{T} = Entry{T}("-", "Temperature shaping factor")
    Te_sep::Entry{T} = Entry{T}("eV", "Separatrix temperature"; default=80.0, check=x -> @assert x > 0.0 "must be: Te_sep > 0.0")
    Te_ped::Entry{T} = Entry{T}("eV", "Pedestal temperature"; check=x -> @assert x > 0.0 "must be: Te_ped > 0.0")
    Te_core::Entry{T} =
        Entry{T}("eV", "Core temperature (NOTE: `Te_core` can be calculated from `ini.equilibrium.presssure_core`)"; check=x -> @assert x > 0.0 "must be: Te_core > 0.0")
    zeff::Entry{T} = Entry{T}("-", "Effective ion charge"; check=x -> @assert x >= 1.0 "must be: zeff > 1.0")
    rot_core::Entry{T} = Entry{T}(IMAS.core_profiles__profiles_1d, :rotation_frequency_tor_sonic)
    ngrid::Entry{Int} = Entry{Int}("-", "Resolution of the core_profiles grid"; default=101, check=x -> @assert x >= 11 "must be: ngrid >= 11")
    bulk::Entry{Symbol} = Entry{Symbol}("-", "Bulk ion species")
    impurity::Entry{Symbol} = Entry{Symbol}("-", "Impurity ion species")
    helium_fraction::Entry{T} = Entry{T}("-", "Helium density / electron density fraction"; check=x -> @assert 0.0 <= x <= 0.5 "must be: 0.0 <= helium_fraction <= 0.5")
    ejima::Entry{T} = Entry{T}("-", "Ejima coefficient"; default=0.4, check=x -> @assert 0.0 <= x < 1.0 "must be: 0.0 <= ejima < 1.0")
    polarized_fuel_fraction::Entry{T} =
        Entry{T}("-", "Spin polarized fuel fraction"; default=0.0, check=x -> @assert 0.0 <= x <= 1.0 "must be: 0.0 <= polarized_fuel_fraction <= 1.0")
    ITB::FUSEparameters__ITB{T} = FUSEparameters__ITB{T}()
end

Base.@kwdef mutable struct FUSEparameters__pf_active{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :pf_active
    technology::Switch{Symbol} = Switch{Symbol}(FusionMaterials.supported_coil_techs(), "-", "PF coils technology")
end

Base.@kwdef mutable struct FUSEparameters__nb_unit{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :nbi
    power_launched::Entry{T} = Entry{T}("W", "Beam power"; check=x -> @assert x >= 0.0 "must be: power_launched >= 0.0")
    beam_energy::Entry{T} = Entry{T}("eV", "Beam energy"; check=x -> @assert x >= 0.0 "must be: beam_energy >= 0.0")
    beam_mass::Entry{T} = Entry{T}("AU", "Beam mass"; default=2.0, check=x -> @assert x >= 1.0 "must be: beam_mass >= 1.0")
    normalized_tangency_radius::Entry{T} = Entry{T}("-", "Tangency radius normalized to major radius "; default=0.6, check=x -> @assert x < 2.0 "must be: beam_mass >= 1.0")
    beam_current_fraction::Entry{Vector{T}}=  Entry{Vector{Vector{T}}} = Entry{Vector{T}}("-", "Beam current fraction", check=x -> @assert sum(x) <= 1.0 )
    current_direction::Entry{Symbol} =  Switch{Symbol}([:co, :counter], "-", "Direction of beam current relative to plasma current",default=false)
    offaxis::Entry{Bool} = Entry{Bool}("-", "Injection neutral beam off axis",default=false)
    template_beam:: Entry{Symbol} =  Switch{Symbol}([:d3d_co,:d3d_counter,:d3d_offaxis, :nstx,:mast_onaxis,:mast_offaxis,:iter_onaxis,:iter_offaxis])
    efficiency_conversion::Entry{T} = Entry{T}(IMAS.nbi__unit___efficiency, :conversion; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_conversion > 0.0")
    efficiency_transmission::Entry{T} = Entry{T}(IMAS.nbi__unit___efficiency, :transmission; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_transmission > 0.0")
end


beamlet.position.r = ini_nbu.r
beamlet.position.z = ini_nbu.z
beamlet.tangency_radius =ini_nbu.tangency_radius
beamlet.angle = ini_nbu.angle
beamlet.beam_current_fraction = ini_nbu.beam_current_fraction 
beamlet.direction = ini_nbu.direction

Base.@kwdef mutable struct FUSEparameters__ec_launcher{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :ec_launcher
    power_launched::Entry{T} = Entry{T}("W", "EC launched power"; check=x -> @assert x >= 0.0 "must be: power_launched >= 0.0")
    efficiency_conversion::Entry{T} = Entry{T}(IMAS.ec_launchers__beam___efficiency, :conversion; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_conversion > 0.0")
    efficiency_transmission::Entry{T} =
        Entry{T}(IMAS.ec_launchers__beam___efficiency, :transmission; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_transmission > 0.0")
end

Base.@kwdef mutable struct FUSEparameters__ic_antenna{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :ic_antenna
    power_launched::Entry{T} = Entry{T}("W", "IC launched power"; check=x -> @assert x >= 0.0 "must be: power_launched >= 0.0")
    efficiency_conversion::Entry{T} = Entry{T}(IMAS.ic_antennas__antenna___efficiency, :conversion; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_conversion > 0.0")
    efficiency_transmission::Entry{T} =
        Entry{T}(IMAS.ic_antennas__antenna___efficiency, :transmission; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_transmission > 0.0")
    efficiency_coupling::Entry{T} = Entry{T}(IMAS.ic_antennas__antenna___efficiency, :coupling; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_coupling > 0.0")
end

Base.@kwdef mutable struct FUSEparameters__lh_antenna{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :lh_antenna
    power_launched::Entry{T} = Entry{T}("W", "LH launched power"; check=x -> @assert x >= 0.0 "must be: power_launched >= 0.0")
    efficiency_conversion::Entry{T} = Entry{T}(IMAS.lh_antennas__antenna___efficiency, :conversion; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_conversion > 0.0")
    efficiency_transmission::Entry{T} =
        Entry{T}(IMAS.lh_antennas__antenna___efficiency, :transmission; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_transmission > 0.0")
    efficiency_coupling::Entry{T} = Entry{T}(IMAS.lh_antennas__antenna___efficiency, :coupling; default=1.0)
end

Base.@kwdef mutable struct FUSEparameters__pellet_launcher{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :pellet_launcher
    frequency::Entry{T} = Entry{T}("Hz", "Frequency of pellets launched"; check=x -> @assert x >= 0.0 "pellet frequency must be >= 0.0")
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
    nose_hfs_fraction::Entry{Float64} = Entry{Float64}("-", "Relative thickness of the TF nose, expressed as a fraction of high-field side TF leg"; default=0.0)
end

Base.@kwdef mutable struct FUSEparameters__oh{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :oh
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
    material::Switch{Symbol} = Switch{Symbol}(FusionMaterials.all_materials(), "-", "Material of the layer")
    coils_inside::Entry{Union{Int,Vector{Int}}} = Entry{Union{Int,Vector{Int}}}(
        "-",
        "List of coils within this layer";
        check=x -> @assert (typeof(x) <: Int && x > 0) || (length(x) > 0 && minimum(x) > 0) "coils_inside must be > 0"
    )
    shape::Switch{BuildLayerShape} = Switch{BuildLayerShape}(layer_shape_options, "-", "Shape of the layer")
    type::Switch{BuildLayerType} = Switch{BuildLayerType}(layer_type_options, "-", "Type of the layer")
    side::Switch{BuildLayerSide} = Switch{BuildLayerSide}(layer_side_options, "-", "Side of the layer")
end

Base.@kwdef mutable struct FUSEparameters__build{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :build
    layers::ParametersVector{FUSEparameters__build_layer{T}} = ParametersVector{FUSEparameters__build_layer{T}}()
    scale_layers_to_R0::Entry{Bool} = Entry{Bool}("-", "Scale layers thicknesses to center plasma equilibrium (R0) in vacuum vessel"; default=true)
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

Base.@kwdef mutable struct FUSEparameters__hcd{T} <: ParametersInit{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :balance_of_plant
    power_scaling_cost_function::Entry{Function} =
        Entry{Function}("-", "EC, IC, LH, NB power optimization cost function, takes dd as input. Eg. dd -> (1.0 - IMAS.tau_e_thermal(dd) / IMAS.tau_e_h98(dd))")
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
    hcd::FUSEparameters__hcd{T}
    build::FUSEparameters__build{T}
    center_stack::FUSEparameters__center_stack{T} #
    tf::FUSEparameters__tf{T}
    oh::FUSEparameters__oh{T}
    bop::FUSEparameters__balance_of_plant{T}
    requirements::FUSEparameters__requirements{T}
end

function ParametersInits{T}() where {T<:Real}
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
        FUSEparameters__hcd{T}(),
        FUSEparameters__build{T}(),
        FUSEparameters__center_stack{T}(),
        FUSEparameters__tf{T}(),
        FUSEparameters__oh{T}(),
        FUSEparameters__balance_of_plant{T}(),
        FUSEparameters__requirements{T}())

    setup_parameters!(ini)

    return ini
end

function ParametersInits(; kw...)
    return ParametersInits{Float64}(; kw...)
end

###############################
# custom dispatches for build #
###############################
"""
    setproperty!(parameters_build::FUSEparameters__build{T}, field::Symbol, layers::AbstractDict{Symbol,<:Real}) where {T<:Real}

Allows users to initialize layers from a dictionary
"""
function Base.setproperty!(parameters_build::FUSEparameters__build{T}, field::Symbol, layers::AbstractDict{Symbol,<:Real}) where {T<:Real}
    @assert field == :layers
    empty!(parameters_build.layers)
    for (k, (name, thickness)) in enumerate(layers)
        layer = FUSEparameters__build_layer{T}()
        push!(parameters_build.layers, layer)

        # name
        layer.name = replace(string(name), "_" => " ")

        # thickness
        layer.thickness = thickness

        # type
        if occursin("OH", uppercase(layer.name)) && occursin("hfs", lowercase(layer.name))
            layer.type = :oh
        elseif occursin("gap ", lowercase(layer.name))
            layer.type = :gap
        elseif lowercase(layer.name) == "plasma"
            layer.type = :plasma
        elseif occursin("OH", uppercase(layer.name))
            layer.type = :oh
        elseif occursin("TF", uppercase(layer.name))
            layer.type = :tf
        elseif occursin("shield", lowercase(layer.name))
            layer.type = :shield
        elseif occursin("blanket", lowercase(layer.name))
            layer.type = :blanket
        elseif occursin("wall", lowercase(layer.name))
            layer.type = :wall
        elseif occursin("vessel", lowercase(layer.name))
            layer.type = :vessel
        elseif occursin("cryostat", lowercase(layer.name))
            layer.type = :cryostat
            layer.shape = :silo
        end

        # side
        if occursin("hfs", lowercase(layer.name))
            layer.side = :hfs
        elseif occursin("lfs", lowercase(layer.name))
            layer.side = :lfs
        else
            if layer.type == _plasma_
                layer.side = :lhfs
            elseif k < length(layers) / 2
                layer.side = :in
            elseif k > length(layers) / 2
                layer.side = :out
            end
        end
    end
end

function Base.setproperty!(parameters_layer::FUSEparameters__build_layer{T}, field::Symbol, val::Symbol) where {T<:Real}
    par = getfield(parameters_layer, field)

    if field == :material && !ismissing(parameters_layer, :type)
        layer_type = parameters_layer.type

        pretty_layer_type = replace("$layer_type", "_" => "")
        allowed_materials = FusionMaterials.supported_material_list(layer_type)

        if val ∉ allowed_materials
            error("$val is not an allowed material for $(pretty_layer_type) layer type. Acceptable materials are $(join(allowed_materials, ", ")).")
        end
    end

    return setproperty!(par, :value, val)
end

"""
    to_index(layers::Vector{FUSEparameters__build_layer{T}}, name::Symbol) where {T<:Real}

Allows accesing parameters layers by their Symbol
"""
function Base.to_index(layers::Vector{FUSEparameters__build_layer{T}}, name::Symbol) where {T<:Real}
    tmp = findfirst(x -> x.name == replace(string(name), "_" => " "), layers)
    if tmp === nothing
        error("Layer `:$name` not found. Valid ini.build.layers are: $([Symbol(replace(layer.name," " => "_")) for layer in layers])")
    end
    return tmp
end

"""
    dict2par!(dct::AbstractDict, par::ParametersVector{<:FUSEparameters__build_layer})

Custom reading from file of FUSEparameters__build_layer
"""
function SimulationParameters.dict2par!(dct::AbstractDict, par::ParametersVector{<:FUSEparameters__build_layer})
    return parent(par).layers = dct
end

"""
    par2ystr(par::ParametersVector{<:FUSEparameters__build_layer}, txt::Vector{String})

Custom writing to file for FUSEparameters__build_layer
"""
function SimulationParameters.par2ystr(par::ParametersVector{<:FUSEparameters__build_layer}, txt::Vector{String}; show_info::Bool=true, skip_defaults::Bool=false)
    for parameter in par
        p = SimulationParameters.path(parameter)
        sp = SimulationParameters.spath(p)
        depth = (count(".", sp) + count("[", sp) - 1) * 2
        pre = " "^depth
        push!(txt, string(pre, replace(getproperty(parameter, :name, "_each_layer_name_and_thickness_"), " " => "_"), ": ", repr(getproperty(parameter, :thickness, 0.0))))
    end
    return txt
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
    equilibrium.tilt = mxh.tilt
    equilibrium.δ = mxh.δ
    equilibrium.ζ = mxh.ζ
    equilibrium.𝚶 = mxh.𝚶
    equilibrium.twist = mxh.twist
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
    boundary_from = ini.equilibrium.boundary_from
    if boundary_from == :ods
        eqt = dd.equilibrium.time_slice[ini.time.simulation_start]
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
            [ini.equilibrium.𝚶, ini.equilibrium.twist],
            [asin(ini.equilibrium.δ), -ini.equilibrium.ζ])
    else
        error("ini.equilibrium.boundary_from must be one of [:scalars, :rz_points, :MXH_params, :ods]")
    end

    if boundary_from == :ods
        # in case of ODS we have all information to generate MXHboundary
        RX = Float64[x_point.r for x_point in eqt.boundary.x_point]
        ZX = Float64[x_point.z for x_point in eqt.boundary.x_point]
        mxhb = MXHboundary(mxh, ini.equilibrium.xpoints in (:upper, :double), ini.equilibrium.xpoints in (:lower, :double), RX, ZX, pr, pz)
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
    load_ods(ini::ParametersAllInits; error_on_missing_coordinates::Bool=true, time_from_ods::Bool=false)

Load ODSs as specified in `ini.ods.filename`

If `time_from_ods==true` then `dd.global_time` and `ini.time.simulation_start` are set from the last available time in the loaded dd0.

If `time_from_ods==false` then `ini.time.simulation_start` takes precedence, and `dd.global_time` is set from that.

NOTE: supports multiple comma-separated filenames

NOTE: ini.general.dd takes priority over ini.general.ods
"""
function load_ods(ini::ParametersAllInits; error_on_missing_coordinates::Bool=true, time_from_ods::Bool=false)
    if !ismissing(ini.general, :dd)
        return ini.general.dd
    end

    dd = load_ods(ini.ods.filename; error_on_missing_coordinates)

    if time_from_ods
        IMAS.last_global_time(dd)
        ini.time.simulation_start = dd.global_time

    else
        dd.global_time = ini.time.simulation_start
    end

    return dd
end

"""
    load_ods(filenames::String; error_on_missing_coordinates::Bool=true)

Load multiple comma-separated filenames into a single dd
"""
function load_ods(filenames::String; error_on_missing_coordinates::Bool=true)
    return load_ods(strip.(split(filenames, ",")); error_on_missing_coordinates)
end

"""
    load_ods(filenames::Vector{<:AbstractString}; error_on_missing_coordinates::Bool=true)

Load multiple ODSs into a single `dd`
"""
function load_ods(filenames::Vector{<:AbstractString}; error_on_missing_coordinates::Bool=true)
    dd = IMAS.dd()
    for filename in filenames
        filename = replace(filename, r"^__FUSE__" => __FUSE__)
        dd1 = IMAS.json2imas(filename; error_on_missing_coordinates)
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

Load the INI parameters from a JSON file with given `filename`
"""
function json2ini(filename::AbstractString, ini::ParametersAllInits=ParametersInits())
    return SimulationParameters.json2par(filename, ini)
end

"""
    ini2yaml(ini::ParametersAllInits, filename::AbstractString; kw...)

Save the INI parameters to a YAML file with given `filename`

`kw` arguments are passed to the YAML.print function
"""
function ini2yaml(ini::ParametersAllInits, filename::AbstractString; kw...)
    return SimulationParameters.par2yaml(ini, filename; kw...)
end

"""
    yaml2ini(filename::AbstractString, ini::ParametersAllInits=ParametersInits())

Load the INI parameters from a YAML file with given `filename`
"""
function yaml2ini(filename::AbstractString, ini::ParametersAllInits=ParametersInits())
    return SimulationParameters.yaml2par(filename, ini)
end

"""
    ini2dict(ini::ParametersAllInits; kw...)

Convert the INI parameters to a dictionary form
"""
function ini2dict(ini::ParametersAllInits; kw...)
    return SimulationParameters.par2dict(ini; kw...)
end

"""
    dict2ini(dict::AbstractDict, ini::ParametersAllInits=ParametersInits())

Convert dict to INI parameters
"""
function dict2ini(dict::AbstractDict, ini::ParametersAllInits=ParametersInits())
    return SimulationParameters.dict2par!(dict, ini)
end

########
# plot #
########
"""
    plot_ini(ini::ParametersAllInits; time0=global_time(ini), cx=false)

Plots ini time dependent time traces including plasma boundary
"""
@recipe function plot_ini(ini::ParametersAllInits; time0=global_time(ini), cx=false)
    id = IMAS.recipe_id_for_help_plot(ini)
    IMAS.assert_type_and_record_argument(id, Float64, "Time to plot"; time0)
    IMAS.assert_type_and_record_argument(id, Bool, "Plot only cross section"; cx)

    # count number of time-dependent parameters
    if !cx
        N = 0
        for par in SimulationParameters.leaves(ini)
            if typeof(par.value) <: Union{TimeData,Function}
                N += 1
            end
        end
        layout := @layout [N + 1]
        w = max(600, Int(ceil(300 * sqrt(N))))
        h = max(400, Int(ceil(200 * sqrt(N))))
        size --> (w, h)
    end

    time_bkp = ini.time.simulation_start
    try
        ini.time.simulation_start = time0

        # plot equilibrium including x-points
        mxhb = MXHboundary(ini)
        wr = wall_radii(mxhb.mxh.R0, mxhb.mxh.minor_radius, ini.build.plasma_gap)
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

        if !cx
            # plot time dependent parameters
            k = 1
            for par in SimulationParameters.leaves(ini)
                if typeof(par.value) <: Union{TimeData,Function}
                    k += 1
                    @series begin
                        label := ""
                        subplot := k
                        time0 := time0
                        par
                    end
                end
            end
        end

    finally
        ini.time.simulation_start = time_bkp
    end
end

using FusionMaterials: FusionMaterials
import OrderedCollections
import SimulationParameters: SwitchOption

import IMAS: BuildLayerType, _plasma_, _gap_, _oh_, _tf_, _shield_, _blanket_, _wall_, _vessel_, _cryostat_, _divertor_
import IMAS: BuildLayerSide, _lfs_, _lhfs_, _hfs_, _in_, _out_
import IMAS: BuildLayerShape, _offset_, _negative_offset_, _convex_hull_, _princeton_D_exact_, _princeton_D_, _princeton_D_scaled_, _rectangle_, _double_ellipse_, _triple_arc_,
    _miller_, _square_miller_, _spline_, _silo_

const layer_shape_options = Dict(Symbol(string(e)[2:end-1]) => SwitchOption(e, string(e)[2:end-1]) for e in instances(IMAS.BuildLayerShape))
const layer_type_options = Dict(Symbol(string(e)[2:end-1]) => SwitchOption(e, string(e)[2:end-1]) for e in instances(IMAS.BuildLayerType))
const layer_side_options = Dict(Symbol(string(e)[2:end-1]) => SwitchOption(e, string(e)[2:end-1]) for e in instances(IMAS.BuildLayerSide))

Base.@kwdef mutable struct FUSEparameters__general{T} <: ParametersInit where {T<:Real}
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

Base.@kwdef mutable struct FUSEparameters__time{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :time
    pulse_shedule_time_basis::Entry{AbstractRange{Float64}} = Entry{AbstractRange{Float64}}("s", "Time basis used to discretize the pulse schedule")
    simulation_start::Entry{Float64} = Entry{Float64}("s", "Time at which the simulation starts"; default=0.0)
end

Base.@kwdef mutable struct FUSEparameters__equilibrium{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :equilibrium
    B0::Entry{T} = Entry{T}(IMAS.equilibrium__vacuum_toroidal_field, :b0)
    R0::Entry{T} = Entry{T}("m", "Geometric genter of the plasma. NOTE: This also scales the radial build layers."; check=x -> @assert x > 0.0 "must be: R0 > 0.0")
    Z0::Entry{T} = Entry{T}("m", "Z offset of the machine midplane"; default=0.0)
    ϵ::Entry{T} = Entry{T}("-", "Plasma inverse aspect ratio (a/R0). NOTE: This also scales the radial build layers."; check=x -> @assert 0.0 < x < 1.0 "must be: 0.0 < ϵ < 1.0")
    κ::Entry{T} = Entry{T}("-", "Plasma elongation. NOTE: If < 1.0 it defines the fraction of maximum controllable elongation estimate.")
    δ::Entry{T} = Entry{T}(IMAS.equilibrium__time_slice___boundary, :triangularity)
    ζ::Entry{T} = Entry{T}(IMAS.equilibrium__time_slice___boundary, :squareness; default=0.0)
    𝚶::Entry{T} = Entry{T}("-", "Ovality of the plasma boundary for up-down asymmetric plasmas"; default=0.0)
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

Base.@kwdef mutable struct FUSEparameters__core_profiles{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :core_profiles
    greenwald_fraction::Entry{T} =
        Entry{T}("-", "Line average electron density expressed as fraction of Greenwald density"; check=x -> @assert x > 0.0 "must be: greenwald_fraction > 0.0")
    greenwald_fraction_ped::Entry{T} =
        Entry{T}("-", "Pedestal electron density expressed as fraction of Greenwald density"; check=x -> @assert x > 0.0 "must be: greenwald_fraction_ped > 0.0")
    ne_ped::Entry{T} = Entry{T}("m^-3", "Pedestal electron density"; check=x -> @assert x > 0.0 "must be: ne_ped > 0.0")
    w_ped::Entry{T} = Entry{T}("-", "Pedestal width expressed in fraction of ψₙ"; default=0.05, check=x -> @assert x > 0.0 "must be: w_ped > 0.0")
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

Base.@kwdef mutable struct FUSEparameters__pf_active{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :pf_active
    n_coils_inside::Entry{Int} = Entry{Int}("-", "Number of PF coils inside of the TF"; check=x -> @assert x >= 0 "must be: n_coils_inside >= 0")
    n_coils_outside::Entry{Int} = Entry{Int}("-", "Number of PF coils outside of the TF"; check=x -> @assert x >= 0 "must be: n_coils_outside >= 0")
    technology::Switch{Symbol} = Switch{Symbol}(supported_coils_techs, "-", "PF coils technology")
end

Base.@kwdef mutable struct FUSEparameters__tf{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :tf
    n_coils::Entry{Int} = Entry{Int}("-", "Number of TF coils"; check=x -> @assert x >= 0 "must be: n_coils >= 0")
    shape::Switch{BuildLayerShape} = Switch{BuildLayerShape}(layer_shape_options, "-", "Shape of the TF coils"; default=:double_ellipse)
    ripple::Entry{T} =
        Entry{T}("-", "Fraction of toroidal field ripple evaluated at the outermost radius of the plasma chamber"; default=0.01, check=x -> @assert x > 0.0 "must be: ripple > 0.0")
    technology::Switch{Symbol} = Switch{Symbol}(supported_coils_techs, "-", "TF coils technology")
end

Base.@kwdef mutable struct FUSEparameters__oh{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :oh
    n_coils::Entry{Int} = Entry{Int}("-", "Number of OH coils"; check=x -> @assert x >= 0 "must be: n_coils >= 0")
    technology::Switch{Symbol} = Switch{Symbol}(supported_coils_techs, "-", "OH coils technology")
end

Base.@kwdef mutable struct FUSEparameters__center_stack{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :center_stack
    bucked::Entry{Bool} = Entry{Bool}("-", "Flag for bucked boundary conditions between TF and OH (and center plug, if present)"; default=false)
    noslip::Entry{Bool} = Entry{Bool}("-", "Flag for no slip conditions between TF and OH (and center plug, if present)"; default=false)
    plug::Entry{Bool} = Entry{Bool}("-", "Flag for center plug"; default=false)
end

Base.@kwdef mutable struct FUSEparameters__nb_unit{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :nbi
    power_launched::Entry{T} = Entry{T}("W", "Beam power"; check=x -> @assert x >= 0.0 "must be: power_launched >= 0.0")
    beam_energy::Entry{T} = Entry{T}("eV", "Beam energy"; check=x -> @assert x >= 0.0 "must be: beam_energy >= 0.0")
    beam_mass::Entry{T} = Entry{T}("AU", "Beam mass"; default=2.0, check=x -> @assert x >= 1.0 "must be: beam_mass >= 1.0")
    toroidal_angle::Entry{T} = Entry{T}("rad", "Toroidal angle of injection"; default=0.0)
    efficiency_conversion::Entry{T} = Entry{T}(IMAS.nbi__unit___efficiency, :conversion; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_conversion > 0.0")
    efficiency_transmission::Entry{T} = Entry{T}(IMAS.nbi__unit___efficiency, :transmission; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_transmission > 0.0")
end

Base.@kwdef mutable struct FUSEparameters__ec_launcher{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :ec_launcher
    power_launched::Entry{T} = Entry{T}("W", "EC launched power"; check=x -> @assert x >= 0.0 "must be: power_launched >= 0.0")
    efficiency_conversion::Entry{T} = Entry{T}(IMAS.ec_launchers__beam___efficiency, :conversion; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_conversion > 0.0")
    efficiency_transmission::Entry{T} =
        Entry{T}(IMAS.ec_launchers__beam___efficiency, :transmission; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_transmission > 0.0")
end

Base.@kwdef mutable struct FUSEparameters__ic_antenna{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :ic_antenna
    power_launched::Entry{T} = Entry{T}("W", "IC launched power")
    efficiency_conversion::Entry{T} = Entry{T}(IMAS.ic_antennas__antenna___efficiency, :conversion; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_conversion > 0.0")
    efficiency_transmission::Entry{T} =
        Entry{T}(IMAS.ic_antennas__antenna___efficiency, :transmission; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_transmission > 0.0")
    efficiency_coupling::Entry{T} = Entry{T}(IMAS.ic_antennas__antenna___efficiency, :coupling; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_coupling > 0.0")
end

Base.@kwdef mutable struct FUSEparameters__lh_antenna{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :lh_antenna
    power_launched::Entry{T} = Entry{T}("W", "LH launched power")
    efficiency_conversion::Entry{T} = Entry{T}(IMAS.lh_antennas__antenna___efficiency, :conversion; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_conversion > 0.0")
    efficiency_transmission::Entry{T} =
        Entry{T}(IMAS.lh_antennas__antenna___efficiency, :transmission; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_transmission > 0.0")
    efficiency_coupling::Entry{T} = Entry{T}(IMAS.lh_antennas__antenna___efficiency, :coupling; default=1.0)
end

Base.@kwdef mutable struct FUSEparameters__build_layer{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :layer
    name::Entry{String} = Entry{String}("-", "Name of the layer")
    thickness::Entry{Float64} =
        Entry{Float64}("-", "Relative thickness of the layer (layers actual thickness is scaled to match plasma R0)"; check=x -> @assert x >= 0.0 "must be: thickness >= 0.0")
    material::Switch{String} = Switch{String}(
        FusionMaterials.available_materials(["blanket_materials", "shield_materials", "structural_materials", "wall_materials", "plasma_material"]),
        "-",
        "Material of the layer"
    )
    shape::Switch{BuildLayerShape} = Switch{BuildLayerShape}(layer_shape_options, "-", "Shape of the layer")
    type::Switch{BuildLayerType} = Switch{BuildLayerType}(layer_type_options, "-", "Type of the layer")
    side::Switch{BuildLayerSide} = Switch{BuildLayerSide}(layer_side_options, "-", "Side of the layer")
end

Base.@kwdef mutable struct FUSEparameters__build{T} <: ParametersInit where {T<:Real}
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

Base.@kwdef mutable struct FUSEparameters__gasc{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :gasc
    filename::Entry{String} = Entry{String}("-", "Output GASC .json file from which data will be loaded")
    case::Entry{Int} = Entry{Int}("-", "Number of the GASC run to load")
end

Base.@kwdef mutable struct FUSEparameters__ods{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :ods
    filename::Entry{String} = Entry{String}("-", "ODS.json file(s) from which equilibrium is loaded. Multiple comma-separated ODSs can be specified.")
end

Base.@kwdef mutable struct FUSEparameters__requirements{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :requirements
    power_electric_net::Entry{T} = Entry{T}(IMAS.requirements, :power_electric_net; check=x -> @assert x > 0.0 "must be: power_electric_net > 0.0")
    flattop_duration::Entry{T} = Entry{T}(IMAS.requirements, :flattop_duration; check=x -> @assert x > 0.0 "must be: flattop_duration > 0.0")
    log10_flattop_duration::Entry{T} =
        Entry{T}("log10(s)", "Log10 value of the duration of the flattop (use Inf for steady-state). Preferred over `flattop_duration` for optimization studies.")
    tritium_breeding_ratio::Entry{T} = Entry{T}(IMAS.requirements, :tritium_breeding_ratio; check=x -> @assert x > 0.0 "must be: tritium_breeding_ratio > 0.0")
    cost::Entry{T} = Entry{T}(IMAS.requirements, :cost; check=x -> @assert x > 0.0 "must be: cost > 0.0")
    ne_peaking::Entry{T} = Entry{T}(IMAS.requirements, :ne_peaking; check=x -> @assert x > 0.0 "must be: ne_peaking > 0.0")
    q_pol_omp::Entry{T} = Entry{T}(IMAS.requirements, :q_pol_omp; check=x -> @assert x > 0.0 "must be: q_pol_omp > 0.0")
    lh_power_threshold_fraction::Entry{T} = Entry{T}(IMAS.requirements, :lh_power_threshold_fraction; check=x -> @assert x > 0.0 "must be: lh_power_threshold_fraction > 0.0")
    h98y2::Entry{T} = Entry{T}(IMAS.requirements, :h98y2; check=x -> @assert x > 0.0 "must be: h98y2 > 0.0")
    hds03::Entry{T} = Entry{T}(IMAS.requirements, :hds03; check=x -> @assert x > 0.0 "must be: hds03 > 0.0")
    βn::Entry{T} = Entry{T}(IMAS.requirements, :βn; check=x -> @assert x > 0.0 "must be: βn > 0.0")
    coil_j_margin::Entry{T} = Entry{T}(IMAS.requirements, :coil_j_margin; check=x -> @assert x > 0.0 "must be: coil_j_margin > 0.0")
    coil_stress_margin::Entry{T} = Entry{T}(IMAS.requirements, :coil_stress_margin; check=x -> @assert x > 0.0 "must be: coil_j_margin > 0.0")
end

mutable struct ParametersInits{T} <: ParametersAllInits where {T<:Real}
    _parent::WeakRef
    _name::Symbol
    general::FUSEparameters__general{T}
    time::FUSEparameters__time{T}
    gasc::FUSEparameters__gasc{T}
    ods::FUSEparameters__ods{T}
    build::FUSEparameters__build{T}
    equilibrium::FUSEparameters__equilibrium{T}
    core_profiles::FUSEparameters__core_profiles{T}
    pf_active::FUSEparameters__pf_active{T}
    tf::FUSEparameters__tf{T}
    oh::FUSEparameters__oh{T}
    center_stack::FUSEparameters__center_stack{T}
    nb_unit::ParametersVector{FUSEparameters__nb_unit{T}}
    ec_launcher::ParametersVector{FUSEparameters__ec_launcher{T}}
    ic_antenna::ParametersVector{FUSEparameters__ic_antenna{T}}
    lh_antenna::ParametersVector{FUSEparameters__lh_antenna{T}}
    requirements::FUSEparameters__requirements{T}
end

function ParametersInits{T}(; n_nb::Int=0, n_ec::Int=0, n_ic::Int=0, n_lh::Int=0, n_layers::Int=0) where {T<:Real}
    ini = ParametersInits{T}(
        WeakRef(nothing),
        :ini,
        FUSEparameters__general{T}(),
        FUSEparameters__time{T}(),
        FUSEparameters__gasc{T}(),
        FUSEparameters__ods{T}(),
        FUSEparameters__build{T}(),
        FUSEparameters__equilibrium{T}(),
        FUSEparameters__core_profiles{T}(),
        FUSEparameters__pf_active{T}(),
        FUSEparameters__tf{T}(),
        FUSEparameters__oh{T}(),
        FUSEparameters__center_stack{T}(),
        ParametersVector{FUSEparameters__nb_unit{T}}(),
        ParametersVector{FUSEparameters__ec_launcher{T}}(),
        ParametersVector{FUSEparameters__ic_antenna{T}}(),
        ParametersVector{FUSEparameters__lh_antenna{T}}(),
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

    for k in 1:n_ic
        push!(ini.ic_antenna, FUSEparameters__ic_antenna{T}())
    end

    for k in 1:n_lh
        push!(ini.lh_antenna, FUSEparameters__lh_antenna{T}())
    end

    setup_parameters!(ini)

    return ini
end

function ParametersInits(args...; kw...)
    return ParametersInits{Float64}(args...; kw...)
end

"""
    ini2json(ini::ParametersAllInits, filename::AbstractString; kw...)

Save the FUSE parameters to a JSON file with give `filename`
`kw` arguments are passed to the JSON.print function
"""
function ini2json(ini::ParametersAllInits, filename::AbstractString; kw...)
    return SimulationParameters.par2json(ini, filename; kw...)
end

function json2ini(filename::AbstractString)
    return SimulationParameters.json2par(filename, ParametersInits())
end

"""
    ini2yaml(ini::ParametersAllInits, filename::AbstractString; kw...)

Save the FUSE parameters to a YAML file with give `filename`
`kw` arguments are passed to the YAML.print function
"""
function ini2yaml(ini::ParametersAllInits, filename::AbstractString; kw...)
    return SimulationParameters.par2yaml(ini, filename; kw...)
end

function yaml2ini(filename::AbstractString)
    return SimulationParameters.yaml2par(filename, ParametersInits())
end

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
    equilibrium.δ = sin(mxh.s[1])
    equilibrium.ζ = -mxh.s[2]
    return equilibrium.𝚶 = mxh.c[1]
end

"""
    IMAS.MXH(equilibrium::FUSEparameters__equilibrium)

return ini.equilibrium boundary expressed in MHX independenty of how the user input it
"""
function IMAS.MXH(ini::ParametersAllInits)
    if ini.general.init_from == :ods
        dd = load_ODSs_from_string(ini.ods.filename)
    else
        dd = IMAS.dd()
    end
    return IMAS.MXH(ini, dd)
end

function IMAS.MXH(ini::ParametersAllInits, dd::IMAS.dd)
    init_from = ini.general.init_from
    if init_from == :ods
        if !ismissing(dd.equilibrium, :time) && length(dd.equilibrium.time) > 0
            dd.global_time = ini.time.simulation_start
            eqt = dd.equilibrium.time_slice[]
            IMAS.flux_surfaces(eqt)
            dd.equilibrium # to avoid GC?
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
            0.0,
            [ini.equilibrium.𝚶, 0.0],
            [asin(ini.equilibrium.δ), -ini.equilibrium.ζ])
    else
        error("ini.equilibrium.boundary_from must be one of [:scalars, :rz_points, :MXH_params, :ods]")
    end

    return mxh
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
    plot_ini(ini::ParametersAllInits)

Plots ini time dependent time traces including plasma boundary
"""
@recipe function plot_ini(ini::ParametersAllInits)
    N = 0
    for par in SimulationParameters.leaves(ini)
        if typeof(par.value) <: Function
            N += 1
        end
    end

    layout := @layout [N + 1]

    mxh = IMAS.MXH(ini)
    @series begin
        label := ""
        subplot := 1
        aspectratio := :equal
        xlim := (0, mxh.R0 * 2)
        mxh
    end

    k = 1
    for par in SimulationParameters.leaves(ini)
        if typeof(par.value) <: Function
            k += 1
            @series begin
                label := ""
                subplot := k
                par
            end
        end
    end
end

"""
    load_ODSs_from_string(filenames::String)

Load multiple comma-separated filenames into a single dd
"""
function load_ODSs_from_string(filenames::String)
    dd = IMAS.dd()
    for filename in split(filenames, ",")
        filename = replace(filename, r"^__FUSE__" => __FUSE__)
        dd1 = IMAS.json2imas(filename)
        merge!(dd, dd1)
    end
    return dd
end
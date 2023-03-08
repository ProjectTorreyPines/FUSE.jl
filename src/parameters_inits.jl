using FusionMaterials: FusionMaterials
import DataStructures

Base.@kwdef mutable struct FUSEparameters__general{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :general
    casename::Entry{String} = Entry(String, "-", "Mnemonic name of the case being run")
    init_from::Switch{Symbol} = Switch(Symbol, [
            :ods => "Load data from ODS saved in .json format (where possible, and fallback on scalars otherwise)",
            :scalars => "Initialize FUSE run from scalar parameters"
        ], "-", "Initialize run from")
end

Base.@kwdef mutable struct FUSEparameters__material{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :material
    wall::Switch{String} = Switch(String, FusionMaterials.available_materials("wall_materials"), "-", "Material used for the wall"; default="Steel, Stainless 316")
    blanket::Switch{String} = Switch(String, FusionMaterials.available_materials("blanket_materials"), "-", "Material used for blanket coils")
    shield::Switch{String} = Switch(String, FusionMaterials.available_materials("shield_materials"), "-", "Material used for the shield")
end

Base.@kwdef mutable struct FUSEparameters__equilibrium{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :equilibrium
    B0::Entry{T} = Entry(T, IMAS.equilibrium__vacuum_toroidal_field, :b0)
    R0::Entry{T} = Entry(T, "m", "Geometric genter of the plasma. NOTE: This also scales the radial build layers.")
    Z0::Entry{T} = Entry(T, "m", "Z offset of the machine midplane"; default=0.0)
    ϵ::Entry{T} = Entry(T, "-", "Plasma inverse aspect ratio. NOTE: This also scales the radial build layers.")
    κ::Entry{T} = Entry(T, IMAS.equilibrium__time_slice___boundary, :elongation)
    δ::Entry{T} = Entry(T, IMAS.equilibrium__time_slice___boundary, :triangularity)
    ζ::Entry{T} = Entry(T, IMAS.equilibrium__time_slice___boundary, :squareness; default=0.0)
    pressure_core::Entry{T} = Entry(T, "Pa", "On axis pressure")
    ip::Entry{T} = Entry(T, IMAS.equilibrium__time_slice___global_quantities, :ip)
    xpoints_number::Entry{Int} = Entry(Int, "-", "Number of x-points")
    ngrid::Entry{Int} = Entry(Int, "-", "Resolution of the equilibrium grid"; default=129)
    field_null_surface::Entry{T} = Entry(T, "-", "ψn value of the field_null_surface. Disable with 0.0"; default=0.5)
    boundary_from::Switch{Symbol} = Switch(Symbol, [:scalars, :MXH_params, :rz_points, :ods], "-", "The starting r, z boundary taken from")
    MXH_params::Entry{Vector{<:T}} = Entry(Vector{<:T}, "-", "Vector of MXH flats")
    rz_points::Entry{Vector{Vector{<:T}}} = Entry(Vector{Vector{<:T}}, "m", "R_Z boundary as Vector{Vector{<:Real}}} : r = rz_points[1], z = rz_points[2]")
end

Base.@kwdef mutable struct FUSEparameters__core_profiles{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :core_profiles
    greenwald_fraction::Entry{T} = Entry(T, "-", "Greenwald fraction, ne_vol / ne_gw")
    ne_ped::Entry{T} = Entry(T, "m^-3", "Pedestal electron density")
    w_ped::Entry{T} = Entry(T, "-", "Pedestal width expressed in fraction of ψₙ", default=0.05)
    T_shaping::Entry{T} = Entry(T, "-", "Temperature shaping factor")
    zeff::Entry{T} = Entry(T, "-", "Effective ion charge")
    rot_core::Entry{T} = Entry(T, IMAS.core_profiles__profiles_1d, :rotation_frequency_tor_sonic)
    ngrid::Entry{Int} = Entry(Int, "-", "Resolution of the core_profiles grid"; default=101)
    bulk::Entry{Symbol} = Entry(Symbol, "-", "Bulk ion species")
    impurity::Entry{Symbol} = Entry(Symbol, "-", "Impurity ion species")
    helium_fraction::Entry{T} = Entry(T, "-", "Helium density / electron density fraction")
    ejima::Entry{T} = Entry(T, "-", "Ejima coefficient")
    polarized_fuel_fraction::Entry{T} = Entry(T, "-", "Spin polarized fuel fraction", default=0.0)
end

Base.@kwdef mutable struct FUSEparameters__coil_tech{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :coil_tech
    material = Switch(String, FusionMaterials.available_materials("magnet_materials"), "-", "Technology used for the coil.")
    temperature::Entry{T} = Entry(T, "K", "Coil temperature")
    thermal_strain::Entry{T} = Entry(T, "-", "Fraction of thermal expansion strain over maximum total strain on coil")
    JxB_strain::Entry{T} = Entry(T, "-", "Fraction of maximum JxB strain over maximum total strain on coil")
    fraction_stainless::Entry{T} = Entry(T, "-", "Fraction of stainless steel in the coil cross-sectional areas")
    fraction_void::Entry{T} = Entry(T, "-", "Fraction of `void` in the coil cross-sectional area. Void is everything (like coolant) that is not structural nor conductor.")
    ratio_SC_to_copper::Entry{T} = Entry(T, "-", "Fraction of superconductor to copper cross-sectional areas")
end

Base.@kwdef mutable struct FUSEparameters__pf_active{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :pf_active
    n_oh_coils::Entry{Int} = Entry(Int, "-", "Number of OH coils")
    n_pf_coils_inside::Entry{Int} = Entry(Int, "-", "Number of PF coils inside of the TF")
    n_pf_coils_outside::Entry{Int} = Entry(Int, "-", "Number of PF coils outside of the TF")
    technology::FUSEparameters__coil_tech{T} = FUSEparameters__coil_tech{T}()
end

Base.@kwdef mutable struct FUSEparameters__tf{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :tf
    n_coils::Entry{Int} = Entry(Int, "-", "Number of TF coils")
    shape::Switch{Symbol} = Switch(Symbol, [:princeton_D_exact, :princeton_D, :princeton_D_scaled, :rectangle, :double_ellipse, :triple_arc, :miller, :spline], "-", "Shape of the TF coils"; default=:princeton_D_scaled)
    ripple::Entry{T} = Entry(T, "-", "Fraction of toroidal field ripple evaluated at the outermost radius of the plasma chamber"; default=0.01)
    technology::FUSEparameters__coil_tech{T} = FUSEparameters__coil_tech{T}()
end


Base.@kwdef mutable struct FUSEparameters__oh{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :oh
    technology::FUSEparameters__coil_tech{T} = FUSEparameters__coil_tech{T}()
end

Base.@kwdef mutable struct FUSEparameters__center_stack{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :center_stack
    bucked::Entry{Bool} = Entry(Bool, "-", "flag for bucked boundary conditions between TF and OH (and center plug, if present)"; default=false)
    noslip::Entry{Bool} = Entry(Bool, "-", "flag for no slip conditions between TF and OH (and center plug, if present)"; default=false)
    plug::Entry{Bool} = Entry(Bool, "-", "flag for center plug"; default=false)
end

Base.@kwdef mutable struct FUSEparameters__nbi{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :nbi
    power_launched::Entry{Union{T,Vector{<:T}}} = Entry(Union{T,Vector{<:T}}, "W", "Beam power")
    beam_energy::Entry{Union{T,Vector{<:T}}} = Entry(Union{T,Vector{<:T}}, "eV", "Beam energy")
    beam_mass::Entry{Union{T,Vector{<:T}}} = Entry(Union{T,Vector{<:T}}, "AU", "Beam mass"; default=2.0)
    toroidal_angle::Entry{Union{T,Vector{<:T}}} = Entry(Union{T,Vector{<:T}}, "rad", "toroidal angle of injection"; default=0.0)
    efficiency_conversion::Entry{Union{T,Vector{<:T}}} = Entry(Union{T,Vector{<:T}}, IMAS.nbi__unit___efficiency, :conversion)
    efficiency_transmission::Entry{Union{T,Vector{<:T}}} = Entry(Union{T,Vector{<:T}}, IMAS.nbi__unit___efficiency, :transmission)
end

Base.@kwdef mutable struct FUSEparameters__ec_launchers{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :ec_launchers
    power_launched::Entry{Union{T,Vector{<:T}}} = Entry(Union{T,Vector{<:T}}, "W", "EC launched power")
    efficiency_conversion::Entry{Union{T,Vector{<:T}}} = Entry(Union{T,Vector{<:T}}, IMAS.ec_launchers__beam___efficiency, :conversion)
    efficiency_transmission::Entry{Union{T,Vector{<:T}}} = Entry(Union{T,Vector{<:T}}, IMAS.ec_launchers__beam___efficiency, :transmission)
end

Base.@kwdef mutable struct FUSEparameters__ic_antennas{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :ic_antennas
    power_launched::Entry{Union{T,Vector{<:T}}} = Entry(Union{T,Vector{<:T}}, "W", "IC launched power")
    efficiency_conversion::Entry{Union{T,Vector{<:T}}} = Entry(Union{T,Vector{<:T}}, IMAS.ic_antennas__antenna___efficiency, :conversion)
    efficiency_transmission::Entry{Union{T,Vector{<:T}}} = Entry(Union{T,Vector{<:T}}, IMAS.ic_antennas__antenna___efficiency, :transmission)
    efficiency_coupling::Entry{Union{T,Vector{<:T}}} = Entry(Union{T,Vector{<:T}}, IMAS.ic_antennas__antenna___efficiency, :coupling)
end

Base.@kwdef mutable struct FUSEparameters__lh_antennas{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :lh_antennas
    power_launched::Entry{Union{T,Vector{<:T}}} = Entry(Union{T,Vector{<:T}}, "W", "LH launched power")
    efficiency_conversion::Entry{Union{T,Vector{<:T}}} = Entry(Union{T,Vector{<:T}}, IMAS.lh_antennas__antenna___efficiency, :conversion)
    efficiency_transmission::Entry{Union{T,Vector{<:T}}} = Entry(Union{T,Vector{<:T}}, IMAS.lh_antennas__antenna___efficiency, :transmission)
    efficiency_coupling::Entry{Union{T,Vector{<:T}}} = Entry(Union{T,Vector{<:T}}, IMAS.lh_antennas__antenna___efficiency, :coupling)
end

Base.@kwdef mutable struct FUSEparameters__build{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :build
    layers::Entry{DataStructures.OrderedDict} = Entry(DataStructures.OrderedDict, "m", "Sorted dictionary of layers thicknesses in radial build")
    blanket::Entry{T} = Entry(T, "-", "Fraction of blanket in radial build")
    shield::Entry{T} = Entry(T, "-", "Fraction of shield in radial build")
    vessel::Entry{T} = Entry(T, "-", "Fraction of vessel in radial build")
    plasma_gap::Entry{T} = Entry(T, "-", "Fraction of vacuum gap between first wall and plasma separatrix in radial build"; default=0.1)
    symmetric::Entry{Bool} = Entry(Bool, "-", "Is the build up-down symmetric")
    n_first_wall_conformal_layers::Entry{Int} = Entry(Int, "-", "Number of layers that are conformal to the first wall"; default=1)
end

Base.@kwdef mutable struct FUSEparameters__gasc{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :gasc
    filename::Entry{String} = Entry(String, "-", "Output GASC .json file from which data will be loaded")
    case::Entry{Int} = Entry(Int, "-", "Number of the GASC run to load")
end

Base.@kwdef mutable struct FUSEparameters__ods{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :ods
    filename::Entry{String} = Entry(String, "-", "ODS.json file from which equilibrium is loaded")
end

Base.@kwdef mutable struct FUSEparameters__requirements{T} <: ParametersInit where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :requirements
    power_electric_net::Entry{T} = Entry(T, "W", "Target net electric power generated by the fusion power plant")
    flattop_duration::Entry{T} = Entry(T, "s", "Target duration of the flattop (use Inf for steady-state)")
    tritium_breeding_ratio::Entry{T} = Entry(T, "-", "Target tritium breeding ratio of the whole plant")
    cost::Entry{T} = Entry(T, "\$M", "Target total FPP cost")
end

mutable struct ParametersInits{T} <: ParametersAllInits where {T<:Real}
    _parent::WeakRef
    _name::Symbol
    general::FUSEparameters__general{T}
    gasc::FUSEparameters__gasc{T}
    ods::FUSEparameters__ods{T}
    build::FUSEparameters__build{T}
    material::FUSEparameters__material{T}
    equilibrium::FUSEparameters__equilibrium{T}
    core_profiles::FUSEparameters__core_profiles{T}
    pf_active::FUSEparameters__pf_active{T}
    tf::FUSEparameters__tf{T}
    oh::FUSEparameters__oh{T}
    center_stack::FUSEparameters__center_stack{T}
    nbi::FUSEparameters__nbi{T}
    ec_launchers::FUSEparameters__ec_launchers{T}
    ic_antennas::FUSEparameters__ic_antennas{T}
    lh_antennas::FUSEparameters__lh_antennas{T}
    requirements::FUSEparameters__requirements{T}
end

function ParametersInits{T}() where {T<:Real}
    ini = ParametersInits{T}(
        WeakRef(nothing),
        :ini,
        FUSEparameters__general{T}(),
        FUSEparameters__gasc{T}(),
        FUSEparameters__ods{T}(),
        FUSEparameters__build{T}(),
        FUSEparameters__material{T}(),
        FUSEparameters__equilibrium{T}(),
        FUSEparameters__core_profiles{T}(),
        FUSEparameters__pf_active{T}(),
        FUSEparameters__tf{T}(),
        FUSEparameters__oh{T}(),
        FUSEparameters__center_stack{T}(),
        FUSEparameters__nbi{T}(),
        FUSEparameters__ec_launchers{T}(),
        FUSEparameters__ic_antennas{T}(),
        FUSEparameters__lh_antennas{T}(),
        FUSEparameters__requirements{T}())
    setup_parameters!(ini)
    return ini
end

function ParametersInits()
    return ParametersInits{Float64}()
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
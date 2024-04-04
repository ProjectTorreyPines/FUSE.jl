import SimulationParameters: SwitchOption

import IMAS: BuildLayerType, _plasma_, _gap_, _oh_, _tf_, _shield_, _blanket_, _wall_, _vessel_, _cryostat_, _divertor_
import IMAS: BuildLayerSide, _lfs_, _lhfs_, _hfs_, _in_, _out_
import IMAS: BuildLayerShape, _offset_, _negative_offset_, _convex_hull_, _princeton_D_exact_, _princeton_D_, _princeton_D_scaled_, _rectangle_, _double_ellipse_,
    _rectangle_ellipse_, _triple_arc_, _miller_, _square_miller_, _spline_, _silo_

const layer_shape_options = Dict(Symbol(string(e)[2:end-1]) => SwitchOption(e, string(e)[2:end-1]) for e in instances(IMAS.BuildLayerShape))
const layer_type_options = Dict(Symbol(string(e)[2:end-1]) => SwitchOption(e, string(e)[2:end-1]) for e in instances(IMAS.BuildLayerType))
const layer_side_options = Dict(Symbol(string(e)[2:end-1]) => SwitchOption(e, string(e)[2:end-1]) for e in instances(IMAS.BuildLayerSide))

Base.@kwdef mutable struct FUSEparameters__tf{T} <: ParametersInitBuild{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :tf
    n_coils::Entry{Int} = Entry{Int}("-", "Number of TF coils"; check=x -> @assert x >= 0 "must be: n_coils >= 0")
    shape::Switch{BuildLayerShape} = Switch{BuildLayerShape}(layer_shape_options, "-", "Shape of the TF coils"; default=:double_ellipse)
    ripple::Entry{T} =
        Entry{T}("-", "Fraction of toroidal field ripple evaluated at the outermost radius of the plasma chamber"; default=0.01, check=x -> @assert x > 0.0 "must be: ripple > 0.0")
    technology::Switch{Symbol} = Switch{Symbol}(FusionMaterials.supported_coil_techs(), "-", "TF coils technology")
end

Base.@kwdef mutable struct FUSEparameters__oh{T} <: ParametersInitBuild{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :oh
    n_coils::Entry{Int} = Entry{Int}("-", "Number of OH coils"; check=x -> @assert x >= 0 "must be: n_coils >= 0")
    technology::Switch{Symbol} = Switch{Symbol}(FusionMaterials.supported_coil_techs(), "-", "OH coils technology")
end

Base.@kwdef mutable struct FUSEparameters__center_stack{T} <: ParametersInitBuild{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :center_stack
    bucked::Entry{Bool} = Entry{Bool}("-", "Flag for bucked boundary conditions between TF and OH (and center plug, if present)"; default=false)
    noslip::Entry{Bool} = Entry{Bool}("-", "Flag for no slip conditions between TF and OH (and center plug, if present)"; default=false)
    plug::Entry{Bool} = Entry{Bool}("-", "Flag for center plug"; default=false)
end

Base.@kwdef mutable struct FUSEparameters__build_layer{T} <: ParametersInitBuild{T}
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

Base.@kwdef mutable struct FUSEparameters__build{T} <: ParametersInitBuild{T}
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

Base.@kwdef mutable struct FUSEparameters__requirements{T} <: ParametersInitBuild{T}
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

mutable struct ParametersInitsBuild{T<:Real} <: ParametersAllInits{T}
    _parent::WeakRef
    _name::Symbol
    general::FUSEparameters__general{T}
    time::FUSEparameters__time{T}
    ods::FUSEparameters__ods{T}
    equilibrium::FUSEparameters__equilibrium{T}
    core_profiles::FUSEparameters__core_profiles{T}
    pf_active::FUSEparameters__pf_active{T}
    nb_unit::ParametersVector{FUSEparameters__nb_unit{T}}
    ec_launcher::ParametersVector{FUSEparameters__ec_launcher{T}}
    pellet_launcher::ParametersVector{FUSEparameters__pellet_launcher{T}}
    ic_antenna::ParametersVector{FUSEparameters__ic_antenna{T}}
    lh_antenna::ParametersVector{FUSEparameters__lh_antenna{T}}
    # below are the machine-build related parameters
    build::FUSEparameters__build{T}
    center_stack::FUSEparameters__center_stack{T} #
    tf::FUSEparameters__tf{T}
    oh::FUSEparameters__oh{T}
    requirements::FUSEparameters__requirements{T}
end

function ParametersInitsBuild{T}(; n_layers::Int=0, kw...) where {T<:Real}

    ini_plasma = ParametersInitsPlasma{T}(; kw...)

    ini = ParametersInitsBuild{T}(
        WeakRef(nothing),
        :ini,
        ini_plasma.general,
        ini_plasma.time,
        ini_plasma.ods,
        ini_plasma.equilibrium,
        ini_plasma.core_profiles,
        ini_plasma.pf_active,
        ini_plasma.nb_unit,
        ini_plasma.ec_launcher,
        ini_plasma.pellet_launcher,
        ini_plasma.ic_antenna,
        ini_plasma.lh_antenna,
        FUSEparameters__build{T}(),
        FUSEparameters__center_stack{T}(),
        FUSEparameters__tf{T}(),
        FUSEparameters__oh{T}(),
        FUSEparameters__requirements{T}())

    for k in 1:n_layers
        push!(ini.build.layers, FUSEparameters__build_layer{T}())
    end

    setup_parameters!(ini)

    return ini
end
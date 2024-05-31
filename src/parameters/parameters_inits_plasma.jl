using FusionMaterials: FusionMaterials

Base.@kwdef mutable struct FUSEparameters__time{T} <: ParametersInitPlasma{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :time
    pulse_shedule_time_basis::Entry{AbstractRange{Float64}} = Entry{AbstractRange{Float64}}("s", "Time basis used to discretize the pulse schedule")
    simulation_start::Entry{Float64} = Entry{Float64}("s", "Time at which the simulation starts"; default=0.0)
end

Base.@kwdef mutable struct FUSEparameters__general{T} <: ParametersInitPlasma{T}
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

Base.@kwdef mutable struct FUSEparameters__equilibrium{T} <: ParametersInitPlasma{T}
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

Base.@kwdef mutable struct FUSEparameters__core_profiles{T} <: ParametersInitPlasma{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :core_profiles
    ne_value::Entry{T} = Entry{T}("-", "Value based on setup method"; check=x -> @assert x > 0.0 "must be > 0.0")
    ne_setting ::Switch{Symbol} = Switch{Symbol}([:ne_ped, :ne_line, :greenwald_fraction, :greenwald_fraction_ped], "-", "Way to set the electron density")
    w_ped::Entry{T} = Entry{T}("-", "Pedestal width expressed in fraction of ψₙ"; default=0.05, check=x -> @assert x > 0.0 "must be: w_ped > 0.0")
    ne_sep_to_ped_ratio::Entry{T} = Entry{T}("-", "Ratio used to set the sepeartrix density based on the pedestal density"; default=0.25, check=x -> @assert x > 0.0 "must be: ne_sep_to_ped_ratio > 0.0")
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

Base.@kwdef mutable struct FUSEparameters__pf_active{T} <: ParametersInitPlasma{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :pf_active
    n_coils_inside::Entry{Int} = Entry{Int}("-", "Number of PF coils inside of the TF"; check=x -> @assert x >= 0 "must be: n_coils_inside >= 0")
    n_coils_outside::Entry{Int} = Entry{Int}("-", "Number of PF coils outside of the TF"; check=x -> @assert x >= 0 "must be: n_coils_outside >= 0")
    technology::Switch{Symbol} = Switch{Symbol}(FusionMaterials.supported_coil_techs(), "-", "PF coils technology")
end

Base.@kwdef mutable struct FUSEparameters__nb_unit{T} <: ParametersInitPlasma{T}
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

Base.@kwdef mutable struct FUSEparameters__ec_launcher{T} <: ParametersInitPlasma{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :ec_launcher
    power_launched::Entry{T} = Entry{T}("W", "EC launched power"; check=x -> @assert x >= 0.0 "must be: power_launched >= 0.0")
    rho_0::Entry{T} = Entry{T}("-", "Desired radial location of the deposition profile"; default=0.5, check=x -> @assert x >= 0.0 "must be: rho_0 >= 0.0")
    width::Entry{T} = Entry{T}("-", "Desired width of the deposition profile"; default=0.025, check=x -> @assert x > 0.0 "must be: width > 0.0")
    efficiency_conversion::Entry{T} = Entry{T}(IMAS.ec_launchers__beam___efficiency, :conversion; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_conversion > 0.0")
    efficiency_transmission::Entry{T} =
        Entry{T}(IMAS.ec_launchers__beam___efficiency, :transmission; default=1.0, check=x -> @assert x > 0.0 "must be: efficiency_transmission > 0.0")
end

Base.@kwdef mutable struct FUSEparameters__ic_antenna{T} <: ParametersInitPlasma{T}
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

Base.@kwdef mutable struct FUSEparameters__lh_antenna{T} <: ParametersInitPlasma{T}
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

Base.@kwdef mutable struct FUSEparameters__pellet_launcher{T} <: ParametersInitPlasma{T}
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

Base.@kwdef mutable struct FUSEparameters__ods{T} <: ParametersInitPlasma{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :ods
    filename::Entry{String} = Entry{String}("-", "ODS.json file(s) from which equilibrium is loaded. Multiple comma-separated ODSs can be specified.")
end

mutable struct ParametersInitsPlasma{T<:Real} <: ParametersAllInits{T}
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
end

function ParametersInitsPlasma{T}(; n_nb::Int=0, n_ec::Int=0, n_pl::Int=0, n_ic::Int=0, n_lh::Int=0) where {T<:Real}
    ini = ParametersInitsPlasma{T}(
        WeakRef(nothing),
        :ini,
        FUSEparameters__general{T}(),
        FUSEparameters__time{T}(),
        FUSEparameters__ods{T}(),
        FUSEparameters__equilibrium{T}(),
        FUSEparameters__core_profiles{T}(),
        FUSEparameters__pf_active{T}(),
        ParametersVector{FUSEparameters__nb_unit{T}}(),
        ParametersVector{FUSEparameters__ec_launcher{T}}(),
        ParametersVector{FUSEparameters__pellet_launcher{T}}(),
        ParametersVector{FUSEparameters__ic_antenna{T}}(),
        ParametersVector{FUSEparameters__lh_antenna{T}}())

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

##############################
# functions operating on ini #
##############################
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
        dd = load_ods(ini)
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
    load_ods(ini::ParametersAllInits)

Load ODSs as specified in `ini.ods.filename`
and sets `dd.global_time` equal to `ini.time.simulation_start`

NOTE: supports multiple comma-separated filenames
"""
function load_ods(ini::ParametersAllInits)
    dd = load_ods(ini.ods.filename)
    dd.global_time = ini.time.simulation_start
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

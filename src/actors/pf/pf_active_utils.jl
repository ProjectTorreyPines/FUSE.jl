import VacuumFields

options_green_model = [
    :point => "one filament per coil",
    :simple => "like :point, but OH coils have three filaments",
    :corners => "like :simple, but PF coils have filaments at the four corners",
    :realistic => "possibly hundreds of filaments per coil (very slow!)"
]

#= ==================================== =#
#  IMAS.pf_active__coil to VacuumFields  #
#= ==================================== =#
mutable struct GS3_IMAS_pf_active__coil{T<:Real,C<:Real} <: VacuumFields.AbstractCoil{T,C}
    imas::IMAS.pf_active__coil{T}
    tech::IMAS.build__pf_active__technology{T}
    time0::Float64
    green_model::Symbol
end

function GS3_IMAS_pf_active__coil(
    pfcoil::IMAS.pf_active__coil{T},
    oh_pf_coil_tech::Union{IMAS.build__oh__technology{T},IMAS.build__pf_active__technology{T}},
    green_model::Symbol) where {T<:Real}

    # convert everything to build__pf_active__technology so that the `coil_tech`
    # type in GS3_IMAS_pf_active__coil is defined at compile time
    coil_tech = IMAS.build__pf_active__technology{T}()
    for field in keys(oh_pf_coil_tech)
        setproperty!(coil_tech, field, getproperty(oh_pf_coil_tech, field))
    end

    return GS3_IMAS_pf_active__coil{T,T}(
        pfcoil,
        coil_tech,
        IMAS.global_time(pfcoil),
        green_model)
end

function IMAS_pf_active__coils(dd::IMAS.dd{D}; green_model::Symbol) where {D<:Real}
    coils = GS3_IMAS_pf_active__coil{D,D}[]
    for coil in dd.pf_active.coil
        if IMAS.is_ohmic_coil(coil)
            coil_tech = dd.build.oh.technology
        else
            coil_tech = dd.build.pf_active.technology
        end
        imas_pf_active__coil = GS3_IMAS_pf_active__coil(coil, coil_tech, green_model)
        push!(coils, imas_pf_active__coil)
    end
    return coils
end

function imas(coil::GS3_IMAS_pf_active__coil)
    return getfield(coil, :imas)
end

function Base.getproperty(coil::GS3_IMAS_pf_active__coil{T}, field::Symbol)::T where {T<:Real}
    pfcoil = getfield(coil, :imas)
    if field ∈ (:r, :z, :width, :height)
        value = getfield(pfcoil.element[1].geometry.rectangle, field)
    elseif field == :turns_with_sign
        value = getfield(pfcoil.element[1], field)
    elseif field == :current
        value = IMAS.get_time_array(pfcoil.current, :data, coil.time0)
    else
        value = getfield(coil, field)
    end
    return value
end

function Base.setproperty!(coil::GS3_IMAS_pf_active__coil, field::Symbol, value::Real)
    pfcoil = getfield(coil, :imas)
    if field == :current
        return IMAS.set_time_array(pfcoil.current, :data, coil.time0, value)
    elseif field ∈ (:r, :z, :width, :height)
        return setfield!(pfcoil.element[1].geometry.rectangle, field, value)
    elseif field == :turns_with_sign
        return setfield!(pfcoil.element[1], field, value)
    else
        setfield!(coil, field, value)
    end
    if field ∈ (:width, :height, :spacing)
        s = sign(getfield(pfcoil, :turns_with_sign))
        turns = Int(ceil(pfcoil.width .* pfcoil.height ./ pfcoil.spacing .^ 2))
        setfield!(pfcoil, :turns_with_sign, s * turns)
    end
    return value
end

#= ==================================== =#
#  IMAS.pf_active__coil to VacuumFields  #
#= ==================================== =#
mutable struct GS_IMAS_pf_active__coil{T<:Real,C<:Real} <: VacuumFields.AbstractCoil{T,C}
    pf_active__coil::IMAS.pf_active__coil{T}
    r::T
    z::T
    width::T
    height::T
    turns_with_sign::T
    spacing::T
    coil_tech::IMAS.build__pf_active__technology{T}
    current_data::Vector{C}
    current_time::Vector{Float64}
    time_index::Int
    green_model::Symbol
end

function GS_IMAS_pf_active__coil(
    pf_active__coil::IMAS.pf_active__coil{T},
    oh_pf_coil_tech::Union{IMAS.build__oh__technology{T},IMAS.build__pf_active__technology{T}},
    green_model::Symbol) where {T<:Real}

    # convert everything to build__pf_active__technology so that the `coil_tech`
    # type in GS_IMAS_pf_active__coil is defined at compile time
    coil_tech = IMAS.build__pf_active__technology{T}()
    for field in keys(oh_pf_coil_tech)
        setproperty!(coil_tech, field, getproperty(oh_pf_coil_tech, field))
    end

    return GS_IMAS_pf_active__coil(
        pf_active__coil,
        pf_active__coil.element[1].geometry.rectangle.r,
        pf_active__coil.element[1].geometry.rectangle.z,
        pf_active__coil.element[1].geometry.rectangle.width,
        pf_active__coil.element[1].geometry.rectangle.height,
        pf_active__coil.element[1].turns_with_sign,
        get_spacing_from_turns(pf_active__coil),
        coil_tech,
        pf_active__coil.current.data,
        pf_active__coil.current.time,
        1,
        green_model)
end

function imas(coil::GS_IMAS_pf_active__coil)
    return getfield(coil, :pf_active__coil)
end

function Base.getproperty(coil::GS_IMAS_pf_active__coil, field::Symbol)
    if field == :current
        tmp = getfield(coil, :current_data)[coil.time_index]
        @assert typeof(tmp) <: Real
        return tmp
    else
        return getfield(coil, field)
    end
end

function Base.setproperty!(coil::GS_IMAS_pf_active__coil, field::Symbol, value::Real)
    if field == :current
        getfield(coil, :current_data)[coil.time_index] = value
    else
        setfield!(coil, field, value)
    end
    if field in (:width, :height, :spacing)
        s = sign(getfield(coil, :turns_with_sign))
        turns = Int(ceil(coil.width .* coil.height ./ coil.spacing .^ 2))
        setfield!(coil, :turns_with_sign, s * turns)
    end
end

function transfer_info_GS_coil_to_IMAS(bd::IMAS.build, coil::GS_IMAS_pf_active__coil)
    pf_active__coil = coil.pf_active__coil
    pf_active__coil.element[1].geometry.rectangle.r = coil.r
    pf_active__coil.element[1].geometry.rectangle.z = coil.z
    pf_active__coil.element[1].geometry.rectangle.width = coil.width
    pf_active__coil.element[1].geometry.rectangle.height = coil.height
    pf_active__coil.element[1].turns_with_sign = float(coil.turns_with_sign)
    pf_active__coil.b_field_max = range(0.1, 30; step=0.1)
    pf_active__coil.temperature = [-1, coil.coil_tech.temperature]
    pf_active__coil.current_limit_max =
        [abs(coil_J_B_crit(b, coil.coil_tech)[1] * area(coil) / coil.turns_with_sign) for b in pf_active__coil.b_field_max, t in pf_active__coil.temperature]
    pf_active__coil.b_field_max_timed.time = coil.current_time
    if IMAS.is_ohmic_coil(pf_active__coil)
        pf_active__coil.b_field_max_timed.data = [bd.oh.max_b_field for time_index in eachindex(coil.current_time)]
    else
        pf_active__coil.b_field_max_timed.data = [coil_selfB(coil, time_index) for time_index in eachindex(coil.current_time)]
    end
    pf_active__coil.current.time = coil.current_time
    return pf_active__coil.current.data = coil.current_data
end

function set_turns_from_spacing!(coil::GS_IMAS_pf_active__coil)
    pf_active__coil = getfield(coil, :pf_active__coil)
    return set_turns_from_spacing!(pf_active__coil, coil.spacing)
end

function set_turns_from_spacing!(pf_active__coil::IMAS.pf_active__coil, spacing::Real)
    s = sign(pf_active__coil.element[1].turns_with_sign)
    return set_turns_from_spacing!(pf_active__coil, spacing, s)
end

function set_turns_from_spacing!(pf_active__coil::IMAS.pf_active__coil, spacing::Real, s::Int)
    return pf_active__coil.element[1].turns_with_sign = float(s * Int(ceil(IMAS.area(pf_active__coil) / spacing^2)))
end

function get_spacing_from_turns(coil::GS_IMAS_pf_active__coil)
    pf_active__coil = getfield(coil, :pf_active__coil)
    return get_spacing_from_turns(pf_active__coil)
end

function get_spacing_from_turns(pf_active__coil::IMAS.pf_active__coil)
    return sqrt((pf_active__coil.element[1].geometry.rectangle.width * pf_active__coil.element[1].geometry.rectangle.height) / abs(pf_active__coil.element[1].turns_with_sign))
end

function area(coil::GS_IMAS_pf_active__coil)
    return IMAS.area(coil.pf_active__coil)
end

"""
    coil_selfB(coil::GS_IMAS_pf_active__coil, time_index::Int=0)

PF coil self-induced magnetic field
NOTE: infinite wire approximation
"""
function coil_selfB(coil::GS_IMAS_pf_active__coil, time_index::Int=0)
    if time_index == 0
        time_index = coil.time_index
    end
    b = abs.(constants.μ_0 * coil.current_data[time_index] * coil.turns_with_sign / (2pi * min(coil.width, coil.height)))
    if b < 0.1
        return 0.1
    else
        return b
    end
end

"""
    fixed_pinned_optim_coils(actor::ActorPFcoilsOpt{D,P}, optimization_scheme::Symbol, coil_struct::Type) where {D<:Real,P<:Real}

Returns tuple of coil_struct (typically `GS_IMAS_pf_active__coil` or `GS3_IMAS_pf_active__coil`) coils organized by their function:

  - fixed: fixed position and current
  - pinned: coils with fixed position but current is optimized
  - optim: coils that have theri position and current optimized
"""
function fixed_pinned_optim_coils(actor::AbstractActor{D,P}, optimization_scheme::Symbol, coil_struct) where {D<:Real,P<:Real}
    dd = actor.dd
    par = actor.par

    fixed_coils = coil_struct{D,D}[]
    pinned_coils = coil_struct{D,D}[]
    optim_coils = coil_struct{D,D}[]
    for coil in dd.pf_active.coil
        if IMAS.is_ohmic_coil(coil)
            coil_tech = dd.build.oh.technology
        else
            coil_tech = dd.build.pf_active.technology
        end
        if coil.identifier == "pinned"
            push!(pinned_coils, coil_struct(coil, coil_tech, par.green_model))
        elseif (coil.identifier == "optim") && (optimization_scheme == :currents)
            push!(pinned_coils, coil_struct(coil, coil_tech, par.green_model))
        elseif coil.identifier == "optim"
            push!(optim_coils, coil_struct(coil, coil_tech, par.green_model))
        elseif coil.identifier == "fixed"
            push!(fixed_coils, coil_struct(coil, coil_tech, par.green_model))
        else
            push!(pinned_coils, coil_struct(coil, coil_tech, par.green_model))
            #     error("$(IMAS.location(coil)).identifier=`$(coil.identifier)` is not valid. Accepted values are [\"optim\", \"pinned\", \"fixed\"]")
        end
    end
    return fixed_coils, pinned_coils, optim_coils
end

#= =================================================================== =#
#  shared between GS_IMAS_pf_active__coil and GS3_IMAS_pf_active__coil  #
#= =================================================================== =#
const All_GS_IMAS_pf_active__coil = Union{GS_IMAS_pf_active__coil,GS3_IMAS_pf_active__coil}

"""
    VacuumFields.Green(coil::All_GS_IMAS_pf_active__coil, R::Real, Z::Real)

Calculates coil green function at given R and Z coordinate
"""
function VacuumFields.Green(coil::All_GS_IMAS_pf_active__coil, R::Real, Z::Real)
    return _gfunc(VacuumFields.Green, coil, R, Z)
end

function VacuumFields.dG_dR(coil::All_GS_IMAS_pf_active__coil, R::Real, Z::Real)
    return _gfunc(VacuumFields.dG_dR, coil, R, Z)
end

function VacuumFields.dG_dZ(coil::All_GS_IMAS_pf_active__coil, R::Real, Z::Real)
    return _gfunc(VacuumFields.dG_dZ, coil, R, Z)
end

function _gfunc(Gfunc::Function, coil::All_GS_IMAS_pf_active__coil, R::Real, Z::Real)
    green_model = getfield(coil, :green_model)
    if green_model == :point # fastest
        return Gfunc(coil.r, coil.z, R, Z, coil.turns_with_sign)

    elseif green_model ∈ (:corners, :simple) # medium
        if IMAS.is_ohmic_coil(imas(coil)) # OH
            n_filaments = max(Int(ceil((coil.height / coil.r) * 2)), 3) # at least 3 filaments, but possibly more as plasma gets closer to the OH
            z_filaments = range(coil.z - (coil.height - coil.width / 2.0) / 2.0, coil.z + (coil.height - coil.width / 2.0) / 2.0; length=n_filaments)
            return sum(Gfunc(coil.r, z, R, Z, coil.turns_with_sign / n_filaments) for z in z_filaments)

        elseif green_model == :simple # PF like point
            return Gfunc(coil.r, coil.z, R, Z, coil.turns_with_sign)

        elseif green_model == :corners # PF with filaments at corners
            return Gfunc(VacuumFields.ParallelogramCoil(coil.r, coil.z, coil.width / 2.0, coil.height / 2.0, 0.0, 90.0, nothing), R, Z, coil.turns_with_sign / 4)

        end

    elseif green_model == :realistic # high-fidelity
        return Gfunc(VacuumFields.ParallelogramCoil(coil.r, coil.z, coil.width, coil.height, 0.0, 90.0, coil.spacing), R, Z)

    else
        error("$(typeof(coil)) green_model can only be (in order of accuracy) :realistic, :corners, :simple, and :point")
    end
end

"""
    encircling_coils(bnd_r::AbstractVector{T}, bnd_z::AbstractVector{T}, r_axis::T, z_axis::T, n_coils::Integer) where {T<:Real}

Generates VacuumFields.PointCoil around the plasma boundary using some educated guesses for where the pf coils should be
"""
function encircling_coils(bnd_r::AbstractVector{T}, bnd_z::AbstractVector{T}, r_axis::T, z_axis::T, n_coils::Integer) where {T<:Real}
    rail_r, rail_z = buffer(bnd_r, bnd_z, (maximum(bnd_r) - minimum(bnd_r)) / 1.5)
    rail_z = (rail_z .- z_axis) .* 1.1 .+ z_axis # give divertors

    valid_r, valid_z = clip_rails(rail_r, rail_z, bnd_r, bnd_z, r_axis, z_axis)

    # normalized distance between -1 and 1
    distance = cumsum(sqrt.(IMAS.gradient(valid_r) .^ 2 .+ IMAS.gradient(valid_z) .^ 2))
    distance = (distance .- distance[1])
    distance = (distance ./ distance[end]) .* 2.0 .- 1.0

    coils_distance = range(-1.0, 1.0, n_coils)
    r_coils = IMAS.interp1d(distance, valid_r).(coils_distance)
    z_coils = IMAS.interp1d(distance, valid_z).(coils_distance)

    r_ohcoils = minimum(bnd_r) / 3.0
    n_oh = Int(ceil((maximum(bnd_z) - minimum(bnd_z)) / r_ohcoils * 2))
    r_ohcoils = 0.01 # this seems to work better

    z_ohcoils = range((minimum(bnd_z) * 2 + minimum(rail_z)) / 3.0, (maximum(bnd_z) * 2 + maximum(rail_z)) / 3.0, n_oh)

    return [
        [VacuumFields.PointCoil(r, z) for (r, z) in zip(z_ohcoils .* 0.0 .+ r_ohcoils, z_ohcoils)];
        [VacuumFields.PointCoil(r, z) for (r, z) in zip(r_coils, z_coils)]
    ]
end

"""
    coil_selfB(coil::IMAS.pf_active__coil{T}, current::T) where {T<:Real}

Evaluates self induced magnetic field of a coil given a current
"""
function coil_selfB(coil::IMAS.pf_active__coil{T}, current::T) where {T<:Real}
    b = abs.(constants.μ_0 * current * coil.element[1].turns_with_sign / (2π * min(coil.element[1].geometry.rectangle.width, coil.element[1].geometry.rectangle.height)))
    if b < 0.1
        return 0.1
    else
        return b
    end
end

"""
    pf_current_limits(pfa::IMAS.pf_active, bd::IMAS.build)

Evaluates the current limit for all PF/OH superconducting coils
"""
function pf_current_limits(pfa::IMAS.pf_active, bd::IMAS.build)
    for coil in pfa.coil
        coil.b_field_max = range(0.1, 30; step=0.1)
        if IMAS.is_ohmic_coil(coil)
            coil_tech = bd.oh.technology
        else
            coil_tech = bd.pf_active.technology
        end
        coil.temperature = [-1, coil_tech.temperature]
        coil.current_limit_max = [abs(coil_J_B_crit(b, coil_tech)[1] * IMAS.area(coil) / coil.element[1].turns_with_sign) for b in coil.b_field_max, t in coil.temperature]
        coil.b_field_max_timed.time = coil.current.time
        if IMAS.is_ohmic_coil(coil)
            if !ismissing(bd.oh, :max_b_field)
                coil.b_field_max_timed.data = [bd.oh.max_b_field for time_index in eachindex(coil.current.time)]
            end
        else
            coil.b_field_max_timed.data = [coil_selfB(coil, coil.current.data[time_index]) for time_index in eachindex(coil.current.time)]
        end
    end
end

@recipe function plot_coil(coil::All_GS_IMAS_pf_active__coil)
    @series begin
        seriestype := :scatter
        marker --> :circle
        label --> ""
        [coil.r], [coil.z]
    end
end

@recipe function plot_coil(coils::AbstractVector{<:All_GS_IMAS_pf_active__coil})
    for (k, coil) in enumerate(coils)
        @series begin
            primary := (k == 1)
            aspect_ratio := :equal
            coil
        end
    end
end

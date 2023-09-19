import VacuumFields

#= ==================================== =#
#  IMAS.pf_active__coil to VacuumFields  #
#= ==================================== =#
mutable struct GS_IMAS_pf_active__coil{T<:Real} <: VacuumFields.AbstractCoil
    pf_active__coil::IMAS.pf_active__coil{T}
    r::T
    z::T
    width::T
    height::T
    turns_with_sign::T
    spacing::T
    coil_tech::IMAS.build__pf_active__technology{T}
    current_data::Vector{T}
    current_time::Vector{T}
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

function Base.getproperty(coil::GS_IMAS_pf_active__coil, field::Symbol)
    if field == :current
        tmp = getfield(coil, :current_data)[coil.time_index]
        @assert typeof(tmp) <: Real
        return tmp
    else
        return getfield(coil, field)
    end
end

function Base.setproperty!(coil::GS_IMAS_pf_active__coil, field::Symbol, value)
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
    if pf_active__coil.name == "OH"
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
    VacuumFields.Green(coil::GS_IMAS_pf_active__coil, R::Real, Z::Real)

Calculates coil green function at given R and Z coordinate
"""
function VacuumFields.Green(coil::GS_IMAS_pf_active__coil, R::Real, Z::Real; n_filaments::Int=3)
    if coil.green_model == :point # fastest
        return VacuumFields.Green(coil.r, coil.z, R, Z, coil.turns_with_sign)

    elseif coil.green_model ∈ (:corners, :simple) # medium
        if coil.pf_active__coil.name == "OH"
            z_filaments = range(coil.z - (coil.height - coil.width / 2.0) / 2.0, coil.z + (coil.height - coil.width / 2.0) / 2.0; length=n_filaments)
            return sum(VacuumFields.Green(coil.r, z, R, Z, coil.turns_with_sign / n_filaments) for z in z_filaments)

        elseif coil.green_model == :corners
            return VacuumFields.Green(VacuumFields.ParallelogramCoil(coil.r, coil.z, coil.width / 2.0, coil.height / 2.0, 0.0, 90.0, nothing), R, Z, coil.turns_with_sign / 4)

        elseif coil.green_model == :simple
            return VacuumFields.Green(coil.r, coil.z, R, Z, coil.turns_with_sign)
        end

    elseif coil.green_model == :realistic # high-fidelity
        return VacuumFields.Green(VacuumFields.ParallelogramCoil(coil.r, coil.z, coil.width, coil.height, 0.0, 90.0, coil.spacing), R, Z)

    else
        error("GS_IMAS_pf_active__coil coil.green_model can only be (in order of accuracy) :realistic, :corners, :simple, and :point")
    end
end

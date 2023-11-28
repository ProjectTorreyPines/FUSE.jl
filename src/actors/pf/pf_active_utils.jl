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
        imas_pf_active__coil.current = 0.0 # initialize coil at global_time
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

"""
    VacuumFields.Green(coil::GS3_IMAS_pf_active__coil, R::Real, Z::Real)

Calculates coil green function at given R and Z coordinate
"""
function VacuumFields.Green(coil::GS3_IMAS_pf_active__coil, R::Real, Z::Real; n_filaments::Int=3)
    return _gfunc(VacuumFields.Green, coil, R, Z; n_filaments)
end
function VacuumFields.dG_dR(coil::GS3_IMAS_pf_active__coil, R::Real, Z::Real; n_filaments::Int=3)
    return _gfunc(VacuumFields.dG_dR, coil, R, Z; n_filaments)
end
function VacuumFields.dG_dZ(coil::GS3_IMAS_pf_active__coil, R::Real, Z::Real; n_filaments::Int=3)
    return _gfunc(VacuumFields.dG_dZ, coil, R, Z; n_filaments)
end

function _gfunc(Gfunc::Function, coil::GS3_IMAS_pf_active__coil, R::Real, Z::Real; n_filaments::Int=3)
    green_model = getfield(coil, :green_model)
    if green_model == :point # fastest
        return Gfunc(coil.r, coil.z, R, Z, coil.turns_with_sign)

    elseif green_model ∈ (:corners, :simple) # medium
        if IMAS.is_ohmic_coil(imas(coil))
            z_filaments = range(coil.z - (coil.height - coil.width / 2.0) / 2.0, coil.z + (coil.height - coil.width / 2.0) / 2.0; length=n_filaments)
            return sum(Gfunc(coil.r, z, R, Z, coil.turns_with_sign / n_filaments) for z in z_filaments)

        elseif green_model == :corners
            return Gfunc(VacuumFields.ParallelogramCoil(coil.r, coil.z, coil.width / 2.0, coil.height / 2.0, 0.0, 90.0, nothing), R, Z, coil.turns_with_sign / 4)

        elseif green_model == :simple
            return Gfunc(coil.r, coil.z, R, Z, coil.turns_with_sign)
        end

    elseif green_model == :realistic # high-fidelity
        return Gfunc(VacuumFields.ParallelogramCoil(coil.r, coil.z, coil.width, coil.height, 0.0, 90.0, coil.spacing), R, Z)

    else
        error("GS3_IMAS_pf_active__coil green_model can only be (in order of accuracy) :realistic, :corners, :simple, and :point")
    end
end

@recipe function plot_coil(C::GS3_IMAS_pf_active__coil)
    @series begin
        seriestype --> :scatter
        marker --> :circle
        markercolor --> :darkgreen
        [C.r], [C.z]
    end
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
    VacuumFields.Green(coil::GS_IMAS_pf_active__coil, R::Real, Z::Real)

Calculates coil green function at given R and Z coordinate
"""
function VacuumFields.Green(coil::GS_IMAS_pf_active__coil, R::Real, Z::Real; n_filaments::Int=3)
    return _gfunc(VacuumFields.Green, coil, R, Z,; n_filaments)
end
function VacuumFields.dG_dR(coil::GS_IMAS_pf_active__coil, R::Real, Z::Real; n_filaments::Int=3)
    return _gfunc(VacuumFields.dG_dR, coil, R, Z,; n_filaments)
end
function VacuumFields.dG_dZ(coil::GS_IMAS_pf_active__coil, R::Real, Z::Real; n_filaments::Int=3)
    return _gfunc(VacuumFields.dG_dZ, coil, R, Z,; n_filaments)
end

function _gfunc(Gfunc::Function, coil::GS_IMAS_pf_active__coil, R::Real, Z::Real; n_filaments::Int=3)
    if coil.green_model == :point # fastest
        return Gfunc(coil.r, coil.z, R, Z, coil.turns_with_sign)

    elseif coil.green_model ∈ (:corners, :simple) # medium
        if IMAS.is_ohmic_coil(coil.pf_active__coil)
            z_filaments = range(coil.z - (coil.height - coil.width / 2.0) / 2.0, coil.z + (coil.height - coil.width / 2.0) / 2.0; length=n_filaments)
            return sum(Gfunc(coil.r, z, R, Z, coil.turns_with_sign / n_filaments) for z in z_filaments)

        elseif coil.green_model == :corners
            return Gfunc(VacuumFields.ParallelogramCoil(coil.r, coil.z, coil.width / 2.0, coil.height / 2.0, 0.0, 90.0, nothing), R, Z, coil.turns_with_sign / 4)

        elseif coil.green_model == :simple
            return Gfunc(coil.r, coil.z, R, Z, coil.turns_with_sign)
        end

    elseif coil.green_model == :realistic # high-fidelity
        return Gfunc(VacuumFields.ParallelogramCoil(coil.r, coil.z, coil.width, coil.height, 0.0, 90.0, coil.spacing), R, Z)

    else
        error("GS_IMAS_pf_active__coil coil.green_model can only be (in order of accuracy) :realistic, :corners, :simple, and :point")
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

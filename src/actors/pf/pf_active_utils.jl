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
mutable struct GS_IMAS_pf_active__coil{T<:Real,C<:Real} <: VacuumFields.AbstractCoil{T,C}
    imas::IMAS.pf_active__coil{T}
    tech::IMAS.build__pf_active__technology{T}
    time0::Float64
    green_model::Symbol
end

function GS_IMAS_pf_active__coil(
    pfcoil::IMAS.pf_active__coil{T},
    oh_pf_coil_tech::Union{IMAS.build__oh__technology{T},IMAS.build__pf_active__technology{T}},
    green_model::Symbol) where {T<:Real}

    # convert everything to build__pf_active__technology so that the `coil_tech`
    # type in GS_IMAS_pf_active__coil is defined at compile time
    coil_tech = IMAS.build__pf_active__technology{T}()
    for field in keys(oh_pf_coil_tech)
        setproperty!(coil_tech, field, getproperty(oh_pf_coil_tech, field))
    end

    return GS_IMAS_pf_active__coil{T,T}(
        pfcoil,
        coil_tech,
        IMAS.global_time(pfcoil),
        green_model)
end

function IMAS_pf_active__coils(dd::IMAS.dd{D}; green_model::Symbol) where {D<:Real}
    coils = GS_IMAS_pf_active__coil{D,D}[]
    for coil in dd.pf_active.coil
        if IMAS.is_ohmic_coil(coil)
            coil_tech = dd.build.oh.technology
        else
            coil_tech = dd.build.pf_active.technology
        end
        imas_pf_active__coil = GS_IMAS_pf_active__coil(coil, coil_tech, green_model)
        push!(coils, imas_pf_active__coil)
    end
    return coils
end

function imas(coil::GS_IMAS_pf_active__coil)
    return getfield(coil, :imas)
end

function Base.getproperty(coil::GS_IMAS_pf_active__coil{T}, field::Symbol) where {T<:Real}
    pfcoil = getfield(coil, :imas)
    if field ∈ (:r, :z, :width, :height)
        value = getfield(pfcoil.element[1].geometry.rectangle, field)
    elseif field == :current
        # IMAS uses current per turn, GS_IMAS_pf_active__coil uses total current
        value = IMAS.get_time_array(pfcoil.current, :data, getfield(coil, :time0)) .* getproperty(pfcoil.element[1], :turns_with_sign, 1.0)
    else
        value = getfield(coil, field)
    end
    return value
end

function Base.setproperty!(coil::GS_IMAS_pf_active__coil, field::Symbol, value::Real)
    pfcoil = getfield(coil, :imas)
    if field == :current
        # IMAS uses current per turn, GS_IMAS_pf_active__coil uses total current
        value = value ./ getproperty(pfcoil.element[1], :turns_with_sign, 1.0)
        return IMAS.set_time_array(pfcoil.current, :data, getfield(coil, :time0), value)
    elseif field ∈ (:r, :z, :width, :height)
        return setfield!(pfcoil.element[1].geometry.rectangle, field, value)
    else
        setfield!(coil, field, value)
    end
    return value
end

"""
    VacuumFields.Green(coil::GS_IMAS_pf_active__coil, R::Real, Z::Real)

Calculates coil green function at given R and Z coordinate
"""
function VacuumFields.Green(coil::GS_IMAS_pf_active__coil, R::Real, Z::Real)
    return _gfunc(VacuumFields.Green, coil, R, Z)
end

function VacuumFields.dG_dR(coil::GS_IMAS_pf_active__coil, R::Real, Z::Real)
    return _gfunc(VacuumFields.dG_dR, coil, R, Z)
end

function VacuumFields.dG_dZ(coil::GS_IMAS_pf_active__coil, R::Real, Z::Real)
    return _gfunc(VacuumFields.dG_dZ, coil, R, Z)
end

function _gfunc(Gfunc::Function, coil::GS_IMAS_pf_active__coil, R::Real, Z::Real)
    green_model = getfield(coil, :green_model)
    if green_model == :point # fastest
        return Gfunc(coil.r, coil.z, R, Z, 1.0)

    elseif green_model ∈ (:corners, :simple) # medium
        if IMAS.is_ohmic_coil(imas(coil)) # OH
            n_filaments = max(Int(ceil((coil.height / coil.r) * 2)), 3) # at least 3 filaments, but possibly more as plasma gets closer to the OH
            z_filaments = range(coil.z - (coil.height - coil.width / 2.0) / 2.0, coil.z + (coil.height - coil.width / 2.0) / 2.0; length=n_filaments)
            return sum(Gfunc(coil.r, z, R, Z, 1.0 / n_filaments) for z in z_filaments)

        elseif green_model == :simple # PF like point
            return Gfunc(coil.r, coil.z, R, Z, 1.0)

        elseif green_model == :corners # PF with filaments at corners
            return Gfunc(VacuumFields.ParallelogramCoil(coil.r, coil.z, coil.width / 2.0, coil.height / 2.0, 0.0, 90.0, nothing), R, Z, 0.25)

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

Evaluates self induced magnetic field of a coil given the total current flowing in it
"""
function coil_selfB(coil::IMAS.pf_active__coil{T}, current::T) where {T<:Real}
    if IMAS.is_ohmic_coil(coil)
        b = IMAS.top_dd(coil).build.oh.max_b_field
    else
        r = sqrt(IMAS.area(coil)) / π
        b = abs.(constants.μ_0 * current / (2π * r))
    end
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
        # OH or PF coil technology
        if IMAS.is_ohmic_coil(coil)
            coil_tech = bd.oh.technology
        else
            coil_tech = bd.pf_active.technology
        end

        # magnetic field of operation
        coil.b_field_max = range(0.1, 30; step=0.1)

        # temperature range of operation
        coil.temperature = [-1, coil_tech.temperature]

        # current limit evaluated at all magnetic fields and temperatures
        coil.current_limit_max = [
            abs(coil_J_B_crit(b, coil_tech)[1] * IMAS.area(coil) * fraction_conductor(coil_tech) / coil.element[1].turns_with_sign) for b in coil.b_field_max,
            t in coil.temperature
        ]

        # maximum magnetic field in time
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

@recipe function plot_coil(coil::GS_IMAS_pf_active__coil)
    @series begin
        seriestype := :scatter
        marker --> :circle
        label --> ""
        [coil.r], [coil.z]
    end
end

@recipe function plot_coil(coils::AbstractVector{<:GS_IMAS_pf_active__coil})
    for (k, coil) in enumerate(coils)
        @series begin
            primary := (k == 1)
            aspect_ratio := :equal
            coil
        end
    end
end

function pack_rail(bd::IMAS.build, λ_regularize::Float64, symmetric::Bool)
    distances = Float64[]
    lbounds = Float64[]
    ubounds = Float64[]
    for rail in bd.pf_active.rail
        if rail.name !== "OH"
            # not symmetric
            if !symmetric
                coil_distances = collect(range(-1.0, 1.0, rail.coils_number + 2))[2:end-1]
                # even symmetric
            elseif mod(rail.coils_number, 2) == 0
                coil_distances = collect(range(-1.0, 1.0, rail.coils_number + 2))[2+Int(rail.coils_number // 2):end-1]
                # odd symmetric
            else
                coil_distances = collect(range(-1.0, 1.0, rail.coils_number + 2))[2+Int((rail.coils_number - 1) // 2)+1:end-1]
            end
            append!(distances, coil_distances)
            append!(lbounds, coil_distances .* 0.0 .- 1.0)
            append!(ubounds, coil_distances .* 0.0 .+ 1.0)
        end
    end
    oh_height_off = Float64[]
    for rail in bd.pf_active.rail
        if rail.name == "OH"
            push!(oh_height_off, 1.0)
            push!(lbounds, 1.0 - 1.0 / rail.coils_number)
            push!(ubounds, 1.0)
            if !symmetric
                push!(oh_height_off, 0.0)
                push!(lbounds, -1.0 / rail.coils_number)
                push!(ubounds, 1.0 / rail.coils_number)
            end
        end
    end
    packed = vcat(distances, oh_height_off, log10(λ_regularize))
    push!(lbounds, -25.0)
    push!(ubounds, -5.0)
    return packed, (lbounds, ubounds)
end

function unpack_rail!(packed::Vector, optim_coils::Vector, symmetric::Bool, bd::IMAS.build)
    λ_regularize = packed[end]
    if symmetric
        n_oh_params = 1
    else
        n_oh_params = 2
    end
    if any(rail.name == "OH" for rail in bd.pf_active.rail)
        oh_height_off = packed[end-n_oh_params:end-1]
        distances = packed[1:end-n_oh_params]
    else
        oh_height_off = Float64[]
        distances = packed[1:end-1]
    end

    if !isempty(optim_coils) # optim_coils have zero length in case of the `static` optimization
        kcoil = 0
        koptim = 0
        koh = 0
        for rail in bd.pf_active.rail
            if rail.name == "OH"
                # mirror OH size when it reaches maximum extent of the rail
                oh_height_off[1] = mirror_bound(oh_height_off[1], 1.0 - 1.0 / rail.coils_number, 1.0)
                if !symmetric
                    offset = mirror_bound(oh_height_off[2], - 1.0 / rail.coils_number, 1.0 / rail.coils_number)
                else
                    offset = 0.0
                end
                z_oh, height_oh = size_oh_coils(rail.outline.z, rail.coils_cleareance, rail.coils_number, oh_height_off[1], offset)
                for k in 1:rail.coils_number
                    koptim += 1
                    koh += 1
                    optim_coils[koptim].z = z_oh[koh]
                    optim_coils[koptim].height = height_oh
                end
            else
                r_interp = IMAS.interp1d(rail.outline.distance, rail.outline.r)
                z_interp = IMAS.interp1d(rail.outline.distance, rail.outline.z)
                # not symmetric
                if !symmetric
                    dkcoil = rail.coils_number
                    coil_distances = distances[kcoil+1:kcoil+dkcoil]
                    # even symmetric
                elseif mod(rail.coils_number, 2) == 0
                    dkcoil = Int(rail.coils_number // 2)
                    coil_distances = distances[kcoil+1:kcoil+dkcoil]
                    coil_distances = vcat(-reverse(coil_distances), coil_distances)
                    # odd symmetric
                else
                    dkcoil = Int((rail.coils_number - 1) // 2)
                    coil_distances = distances[kcoil+1:kcoil+dkcoil]
                    coil_distances = vcat(-reverse(coil_distances), 0.0, coil_distances)
                end
                kcoil += dkcoil

                # mirror coil position when they reach the end of the rail
                coil_distances = mirror_bound.(coil_distances, -1.0, 1.0)

                # get coils r and z from distances
                r_coils = r_interp.(coil_distances)
                z_coils = z_interp.(coil_distances)

                # assign to optim coils
                for k in 1:rail.coils_number
                    koptim += 1
                    optim_coils[koptim].r = r_coils[k]
                    optim_coils[koptim].z = z_coils[k]
                end
            end
        end
    end

    return 10^λ_regularize
end

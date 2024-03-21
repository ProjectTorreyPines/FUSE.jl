import VacuumFields

options_green_model = [
    :point => "One filament per coil",
    :quad => "Quadrilateral coil with quadrature integration"
]

#= ==================================== =#
#  IMAS.pf_active__coil to VacuumFields  #
#= ==================================== =#
mutable struct GS_IMAS_pf_active__coil{T1<:Real,T2<:Real,T3<:Real} <: VacuumFields.AbstractCoil{T1,T2,T3}
    imas::IMAS.pf_active__coil{T1}
    tech::IMAS.build__pf_active__technology{T1}
    time0::Float64
    green_model::Symbol
end

function GS_IMAS_pf_active__coil(
    pfcoil::IMAS.pf_active__coil{T},
    oh_pf_coil_tech::Union{IMAS.build__oh__technology{T},IMAS.build__pf_active__technology{T}},
    green_model::Symbol,
    default_resistance::Float64=1e-6) where {T<:Real}

    # convert everything to build__pf_active__technology so that the `coil_tech`
    # type in GS_IMAS_pf_active__coil is defined at compile time
    coil_tech = IMAS.build__pf_active__technology{T}()
    for field in keys(oh_pf_coil_tech)
        setproperty!(coil_tech, field, getproperty(oh_pf_coil_tech, field))
    end

    coil = GS_IMAS_pf_active__coil{T,T,T}(
            pfcoil,
            coil_tech,
            global_time(pfcoil),
            green_model)

    mat_pf = Material(coil_tech)
    sigma = mat_pf.electrical_conductivity
    if ismissing(mat_pf) || ismissing(sigma)
        coil.resistance = default_resistance
    else
        coil.resistance = VacuumFields.resistance(coil, 1.0 / sigma(temperature=0.0), :parallel)
    end

    return coil

end

function IMAS_pf_active__coils(dd::IMAS.dd{D}; green_model::Symbol=:quad, zero_currents::Bool=false) where {D<:Real}
    coils = GS_IMAS_pf_active__coil{D,D}[]
    for coil in dd.pf_active.coil
        if zero_currents
            @ddtime(coil.current.data = 0.0)   # zero currents for all coils
        end
        if :shaping ∉ [IMAS.index_2_name(coil.function)[f.index] for f in coil.function]
            continue
        end
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
    elseif field == :resistance
        value = getfield(pfcoil, field)
    elseif field == :turns
        value = Int(abs(getproperty(pfcoil.element[1], :turns_with_sign, 1.0)))
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
    elseif field == :resistance
        setfield!(pfcoil, field, value)
    elseif field == :turns
        val = abs(value) * sign(getproperty(pfcoil.element[1], :turns_with_sign, 1.0))
        setfield!(pfcoil.element[1], :turns_with_sign, val)
    else
        setfield!(coil, field, value)
    end
    return value
end

"""
    VacuumFields.Green(coil::GS_IMAS_pf_active__coil, R::Real, Z::Real)

Calculates coil green function at given R and Z coordinate
"""
function VacuumFields.Green(coil::GS_IMAS_pf_active__coil, R::Real, Z::Real, scale_factor::Real=1.0; kwargs...)
    return _gfunc(VacuumFields.Green, coil, R, Z)
end

function VacuumFields.dG_dR(coil::GS_IMAS_pf_active__coil, R::Real, Z::Real, scale_factor::Real=1.0; kwargs...)
    return _gfunc(VacuumFields.dG_dR, coil, R, Z)
end

function VacuumFields.dG_dZ(coil::GS_IMAS_pf_active__coil, R::Real, Z::Real, scale_factor::Real=1.0; kwargs...)
    return _gfunc(VacuumFields.dG_dZ, coil, R, Z)
end

function _gfunc(Gfunc::Function, coil::GS_IMAS_pf_active__coil, R::Real, Z::Real, scale_factor::Real=1.0; xorder::Int=3, yorder::Int=3)
    green_model = getfield(coil, :green_model)

    if green_model == :point # low-fidelity
        oute = IMAS.outline(coil.imas.element[1])
        rc0, zc0 = IMAS.centroid(oute.r, oute.z)
        return Gfunc(rc0, zc0, R, Z)

    elseif green_model == :quad # high-fidelity
        return Gfunc(coil.imas, R, Z)

    else
        error("$(typeof(coil)) green_model is `$(green_model)` but it can only be `:point` or `:quad`")

    end
end


# Mutual inductance

function VacuumFields.mutual(C1::GS_IMAS_pf_active__coil, C2::GS_IMAS_pf_active__coil; xorder::Int=3, yorder::Int=3)

    gm1 = getfield(C1, :green_model)
    gm2 = getfield(C2, :green_model)
    @assert gm1 in (:point, :quad)
    @assert gm2 in (:point, :quad)

    if gm1 === :quad && gm2 === :quad
        VacuumFields.mutual(C1.imas, C2.imas; xorder, yorder)
    else
        fac = -2π * VacuumFields.μ₀ * C1.turns * C2.turns
        if gm1 === :point && gm2 === :point
            return fac * Green(C1.r, C1.z, C2.r, C2.z)
        elseif gm1 === :point
            return fac * Green(C2.imas, C1.r, C1.r)
        else
            return fac * Green(C1.imas, C2.r, C2.r)
        end
    end
end

function VacuumFields.mutual(C1::GS_IMAS_pf_active__coil, C2::VacuumFields.AbstractCoil; xorder::Int=3, yorder::Int=3)

    green_model = getfield(C1, :green_model)
    if green_model == :point # fastest
        fac = -2π * VacuumFields.μ₀ * VacuumFields.turns(C2) * C1.turns
        return fac * Green(C2, C1.r, C1.z)

    elseif green_model == :quad # high-fidelity
        return VacuumFields.mutual(C1.imas, C2; xorder, yorder)

    else
        error("$(typeof(C2)) green_model can only be (in order of accuracy) :quad and :point")
    end
end

function VacuumFields.mutual(C1::VacuumFields.AbstractCoil, C2::GS_IMAS_pf_active__coil; xorder::Int=3, yorder::Int=3)

    green_model = getfield(C2, :green_model)
    if green_model == :point # fastest
        fac = -2π * VacuumFields.μ₀ * VacuumFields.turns(C1) * C2.turns
        return fac * Green(C1, C2.r, C2.z)

    elseif green_model == :quad # high-fidelity
        return VacuumFields.mutual(C1, C2.imas; xorder, yorder)

    else
        error("$(typeof(C2)) green_model can only be (in order of accuracy) :quad and :point")
    end
end

function VacuumFields.mutual(C1::VacuumFields.AbstractCoil, C2::GS_IMAS_pf_active__coil; xorder::Int=3, yorder::Int=3)

    green_model = getfield(C2, :green_model)
    if green_model == :point # fastest
        fac = -2π * VacuumFields.μ₀ * VacuumFields.turns(C1) * C2.turns
        return fac * Green(C1, C2.r, C2.z)

    elseif green_model == :quad # high-fidelity
        return VacuumFields.mutual(C1, C2.imas; xorder, yorder)

    else
        error("$(typeof(C2)) green_model can only be (in order of accuracy) :quad and :point")
    end
end


function VacuumFields._pfunc(Pfunc, image::VacuumFields.Image, C::GS_IMAS_pf_active__coil, δZ;
                COCOS::MXHEquilibrium.COCOS=MXHEquilibrium.cocos(11),
                xorder::Int=3, yorder::Int=3)

    green_model = getfield(C, :green_model)
    if green_model == :point # fastest
        PC = VacuumFields.PointCoil(C.r, C.z; C.turns)
        return VacuumFields._pfunc(Pfunc, image, PC, δZ; COCOS)

    elseif green_model == :quad # high-fidelity
        return VacuumFields._pfunc(Pfunc, image, C.imas, δZ; COCOS, xorder, yorder)

    else
        error("$(typeof(C)) green_model can only be (in order of accuracy) :quad and :point")
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
        [VacuumFields.PointCoil(r, z) for (r, z) in zip(z_ohcoils .* 0.0 .+ r_ohcoils, z_ohcoils)]
        [VacuumFields.PointCoil(r, z) for (r, z) in zip(r_coils, z_coils)]
    ]
end

"""
    coil_selfB(coil::IMAS.pf_active__coil{T}, total_current::T) where {T<:Real}

Evaluates self induced magnetic field of a coil given the TOTAL current flowing in it
"""
function coil_selfB(coil::IMAS.pf_active__coil{T}, total_current::T) where {T<:Real}
    if IMAS.is_ohmic_coil(coil)
        # for the ohmic solenoid we take the maximum b-field used for the flux-swing
        b = IMAS.top_dd(coil).build.oh.max_b_field
    else
        r = sqrt(IMAS.area(coil)) / π
        b = abs.(constants.μ_0 * total_current / (2π * r))
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
        mat_pf = Material(coil_tech)

        # magnetic field of operation
        coil.b_field_max = range(0.1, 30; step=0.1)

        # temperature range of operation
        coil.temperature = [-1, coil_tech.temperature]

        # current limit evaluated at all magnetic fields and temperatures
        coil.current_limit_max = [
            abs(mat_pf.critical_current_density(; Bext=b) * IMAS.area(coil) * IMAS.fraction_conductor(coil_tech) / coil.element[1].turns_with_sign) for
            b in coil.b_field_max,
            t in coil.temperature
        ]

        # maximum magnetic field in time
        coil.b_field_max_timed.time = coil.current.time
        if IMAS.is_ohmic_coil(coil)
            if !ismissing(bd.oh, :max_b_field)
                coil.b_field_max_timed.data = [bd.oh.max_b_field for time_index in eachindex(coil.current.time)]
            end
        else
            coil.b_field_max_timed.data = [coil_selfB(coil, coil.current.data[time_index] .* coil.element[1].turns_with_sign) for time_index in eachindex(coil.current.time)]
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
                coil_distances = collect(range(-1.0, 1.0, rail.coils_number))[1:end]
                # even symmetric
            elseif mod(rail.coils_number, 2) == 0
                coil_distances = collect(range(-1.0, 1.0, rail.coils_number))[1+Int(rail.coils_number // 2):end]
                # odd symmetric
            else
                coil_distances = collect(range(-1.0, 1.0, rail.coils_number))[1+Int((rail.coils_number - 1) // 2)+1:end]
            end
            append!(distances, coil_distances)
            if !symmetric
                append!(lbounds, coil_distances .* 0.0 .- 1.0)
            else
                append!(lbounds, coil_distances .* 0.0)
            end
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
                push!(lbounds, -2.0 / rail.coils_number)
                push!(ubounds, 2.0 / rail.coils_number)
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
                    offset = mirror_bound(oh_height_off[2], -2.0 / rail.coils_number, 2.0 / rail.coils_number)
                else
                    offset = 0.0
                end
                z_oh, height_oh = size_oh_coils(rail.outline.z, rail.coils_cleareance, rail.coils_number, oh_height_off[1], 0.0)
                z_oh = z_oh .+ offset # allow offset to move the whole CS stack independently of the CS rail
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

function size_pf_active(coils::AbstractVector{<:GS_IMAS_pf_active__coil}; tolerance::Float64=0.4, min_size::Float64=1.0)
    function optimal_area(x; coil)
        area = abs(x[1])

        pfcoil = getfield(coil, :imas)

        height = width = sqrt(area)
        pfcoil.element[1].geometry.rectangle.height = height
        pfcoil.element[1].geometry.rectangle.width = width

        mat = Material(coil.tech)
        Bext = coil_selfB(pfcoil, coil.current)

        needed_conductor_area = abs(coil.current) / mat.critical_current_density(; Bext)
        needed_area = needed_conductor_area / IMAS.fraction_conductor(coil.tech) * (1.0 .+ tolerance)

        cost = (area - needed_area)^2
        return cost
    end

    # find optimal area for each coil
    areas = Float64[]
    for coil in coils
        pfcoil = getfield(coil, :imas)
        if !IMAS.is_ohmic_coil(pfcoil)
            res = Optim.optimize(x -> optimal_area(x; coil), [0.1], Optim.NelderMead())
            push!(areas, abs(res.minimizer[1]))
        end
    end

    # set the area of the coils, with a minimum size given by the norm
    msa = norm(areas) / length(areas)
    k = 0
    for coil in coils
        pfcoil = getfield(coil, :imas)
        if !IMAS.is_ohmic_coil(pfcoil)
            k += 1
            optimal_area(max(areas[k], min_size * msa); coil)
        end
    end
end

#= ============================================= =#
#  Visualization of IMAS.pf_active.coil as table  #
#= ============================================= =#
function DataFrames.DataFrame(coils::IMAS.IDSvector{<:IMAS.pf_active__coil})

    df = DataFrames.DataFrame(;
        name=String[],
        var"function"=Vector{Symbol}[],
        n_elements=Int[],
        n_total_turns=Float64[]
    )

    for coil in coils
        func = [IMAS.index_2_name(coil.function)[f.index] for f in coil.function]
        turns = sum(getproperty(element, :turns_with_sign, 1.0) for element in coil.element)
        push!(df, [coil.name, func, length(coil.element), sum(turns)])
    end

    return df
end

function Base.show(io::IO, ::MIME"text/plain", coils::IMAS.IDSvector{<:IMAS.pf_active__coil})
    old_lines = get(ENV, "LINES", missing)
    old_columns = get(ENV, "COLUMNS", missing)
    df = DataFrames.DataFrame(coils)
    try
        ENV["LINES"] = 1000
        ENV["COLUMNS"] = 1000
        return show(io::IO, df)
    finally
        if old_lines === missing
            delete!(ENV, "LINES")
        else
            ENV["LINES"] = old_lines
        end
        if old_columns === missing
            delete!(ENV, "COLUMNS")
        else
            ENV["COLUMNS"] = old_columns
        end
    end
end
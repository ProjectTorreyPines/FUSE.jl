import VacuumFields
import DataFrames

const options_green_model = [
    :point => "One filament per coil",
    :quad => "Quadrilateral coil with quadrature integration"
]

"""
    encircling_coils(bnd_r::AbstractVector{T1}, bnd_z::AbstractVector{T1}, r_axis::T2, z_axis::T2, n_coils::Integer) where {T1<:Real,T2<:Real}

Generates pf_active.coil around the plasma boundary making some educated guesses for where the PF and OH coils should be.

`n_coils` is used to set equal number of "OH" and "PF" coils
"""
function encircling_coils(bnd_r::AbstractVector{T1}, bnd_z::AbstractVector{T1}, r_axis::T2, z_axis::T2, n_coils::Integer) where {T1<:Real,T2<:Real}
    # define PF rail
    a = (maximum(bnd_r) - minimum(bnd_r)) / 2
    rail_r, rail_z = buffer(bnd_r, bnd_z, a * 1.3)
    valid_r, valid_z = clip_rails(rail_r, rail_z, bnd_r, bnd_z, r_axis, z_axis, _lfs_)

    # OH positions and sizes
    r_ohcoils = minimum(bnd_r) / 3
    z_oh_low = (min(valid_z[1], valid_z[end]) + minimum(bnd_z)) / 2.0
    z_oh_high = (max(valid_z[1], valid_z[end]) + maximum(bnd_z)) / 2.0
    z_ohcoils, h_oh = size_oh_coils(z_oh_low, z_oh_high, 0.0, n_coils)
    w_oh = minimum(bnd_r) / 3

    # PF posistions and sizes
    r_coils, z_coils = IMAS.resample_2d_path(valid_r, valid_z; n_points=n_coils, method=:cubic)
    w_pf = h_pf = IMAS.perimeter(valid_r, valid_z) / n_coils / sqrt(2.0)

    # OH coils
    coils = IMAS.pf_active__coil{T1}[]
    for (kk, (r, z)) in enumerate(zip(z_ohcoils .* 0.0 .+ r_ohcoils, z_ohcoils))
        coil = IMAS.pf_active__coil{T1}()
        coil.identifier = "optim"
        coil.name = "OH $kk"
        resize!(coil.element, 1)
        coil.element[1].turns_with_sign = 1.0
        pf_geo = coil.element[1].geometry
        pf_geo.geometry_type = 2
        pf_geo.rectangle.r = r
        pf_geo.rectangle.z = z
        pf_geo.rectangle.width = w_oh
        pf_geo.rectangle.height = h_oh
        func = resize!(coil.function, :shaping)
        func.description = "OH"
        push!(coils, coil)
    end

    # PF coils
    for (kk, (r, z)) in enumerate(zip(r_coils, z_coils))
        coil = IMAS.pf_active__coil{T1}()
        coil.identifier = "optim"
        coil.name = "PF $kk"
        resize!(coil.element, 1)
        coil.element[1].turns_with_sign = 1.0
        pf_geo = coil.element[1].geometry
        pf_geo.geometry_type = 2
        pf_geo.rectangle.r = r
        pf_geo.rectangle.z = z
        pf_geo.rectangle.width = w_pf
        pf_geo.rectangle.height = h_pf
        func = resize!(coil.function, :shaping)
        func.description = "PF"
        push!(coils, coil)
    end

    return coils
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
        b = abs.(IMAS.mks.μ_0 * total_current / (2π * r))
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

        if !ismissing(coil_tech, :material)
            # magnetic field of operation
            coil.b_field_max = range(0.1, 30; step=0.1)

            # temperature range of operation
            coil.temperature = [-1, coil_tech.temperature]

            # current limit evaluated at all magnetic fields and temperatures
            mat_pf = Material(coil_tech)
            coil_area = IMAS.area(coil)
            frac_conductor = IMAS.fraction_conductor(coil_tech)
            turns = coil.element[1].turns_with_sign
            coil.current_limit_max = [
                abs(mat_pf.critical_current_density(; Bext=b) * coil_area * frac_conductor / turns) for
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
end

function pack_rail(bd::IMAS.build, λ_regularize::Float64, symmetric::Bool)
    distances = Float64[]
    lbounds = Float64[]
    ubounds = Float64[]
    for rail in bd.pf_active.rail
        if rail.name == "PF"
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
    if any(rail.name == "OH" for rail in bd.pf_active.rail)
        if symmetric
            n_oh_params = 1
        else
            n_oh_params = 2
        end
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
                oh_height_off[1] = mirror_bound(oh_height_off[1], 0.8, 1.0)
                # allow ± one coil offset
                if symmetric
                    offset = 0.0
                else
                    offset = mirror_bound(oh_height_off[2], -1.0, 1.0)
                end
                z_oh, height_oh = size_oh_coils(minimum(rail.outline.z), maximum(rail.outline.z), rail.coils_cleareance, rail.coils_number, oh_height_off[1], offset)
                r_interp = IMAS.interp1d(rail.outline.distance, rail.outline.r)
                for k in 1:rail.coils_number
                    koptim += 1
                    koh += 1
                    coil_distance = (z_oh[koh] - minimum(rail.outline.z)) / (maximum(rail.outline.z) - minimum(rail.outline.z)) * 2 - 1
                    optim_coils[koptim].z = r_interp(coil_distance)
                    optim_coils[koptim].z = z_oh[koh]
                    optim_coils[koptim].height = height_oh
                end
            elseif rail.name == "PF"
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

function size_pf_active(
    coils::AbstractVector{<:VacuumFields.GS_IMAS_pf_active__coil},
    eqt::IMAS.equilibrium__time_slice;
    tolerance::Float64=0.0,
    min_size::Float64=0.1,
    symmetric::Bool
)
    Rcenter = eqt.boundary.geometric_axis.r
    Zcenter = eqt.boundary.geometric_axis.z

    function optimal_area(x; coil, pfcoil, r0, z0, width0, height0)
        area = abs(x[1])

        height = sqrt(area)
        width = area / height
        pfcoil.element[1].geometry.rectangle.height = height
        pfcoil.element[1].geometry.rectangle.width = width
        pfcoil.element[1].geometry.rectangle.r = r0 + sign(r0 - Rcenter) * (width - width0) / 2.0
        pfcoil.element[1].geometry.rectangle.z = z0 + sign(z0 - Zcenter) * (height - height0) / 2.0

        mat = Material(coil.tech)
        Bext = coil_selfB(pfcoil, coil.current)

        needed_conductor_area = abs(coil.current) / mat.critical_current_density(; Bext)
        needed_area = needed_conductor_area / IMAS.fraction_conductor(coil.tech) * (1.0 .+ tolerance)
        if needed_area > 1E6 # to handle cases where needed_area == Inf
            needed_area = 1E6
        end
        cost = (area - needed_area)^2
        return cost
    end

    # find optimal area for each coil
    areas = Float64[]
    for coil in coils
        pfcoil = getfield(coil, :imas)
        if !IMAS.is_ohmic_coil(pfcoil)
            r0 = pfcoil.element[1].geometry.rectangle.r
            z0 = pfcoil.element[1].geometry.rectangle.z
            width0 = pfcoil.element[1].geometry.rectangle.width
            height0 = pfcoil.element[1].geometry.rectangle.height
            res = Optim.optimize(x -> optimal_area(x; coil, pfcoil, r0, z0, width0, height0), [0.1], Optim.NelderMead())
            push!(areas, abs(res.minimizer[1]))
        end
    end

    # find symmetric coils and make them equal area
    if symmetric
        symmetric_pairs = Tuple{Int,Int}[]
        k1 = 0
        for coil1 in coils
            pfcoil1 = getfield(coil1, :imas)
            if !IMAS.is_ohmic_coil(pfcoil1)
                k1 += 1
                k2 = 0
                min_symmetric_distance = Inf
                for coil2 in coils
                    pfcoil2 = getfield(coil2, :imas)
                    if !IMAS.is_ohmic_coil(pfcoil2)
                        k2 += 1
                        if k1 > k2
                            symmetric_distance = sqrt(
                                (pfcoil1.element[1].geometry.rectangle.r - pfcoil2.element[1].geometry.rectangle.r)^2 +
                                (pfcoil1.element[1].geometry.rectangle.z + pfcoil2.element[1].geometry.rectangle.z)^2
                            )
                            if symmetric_distance < min_symmetric_distance && symmetric_distance < sqrt(max(areas[k2], areas[k1]))
                                push!(symmetric_pairs, (k1, k2))
                                min_symmetric_distance = symmetric_distance
                            end
                        end
                    end
                end
            end
            if !isempty(symmetric_pairs)
                k1, k2 = pop!(symmetric_pairs)
                max_area = max(areas[k2], areas[k1])
                areas[k2] = max_area
                areas[k1] = max_area
            end
        end
    end

    # set the area of the coils, with a minimum size given by the norm
    msa = norm(areas) / length(areas)
    k = 0
    for coil in coils
        pfcoil = getfield(coil, :imas)
        if !IMAS.is_ohmic_coil(pfcoil)
            k += 1
            r0 = pfcoil.element[1].geometry.rectangle.r
            z0 = pfcoil.element[1].geometry.rectangle.z
            width0 = pfcoil.element[1].geometry.rectangle.width
            height0 = pfcoil.element[1].geometry.rectangle.height
            optimal_area(max(areas[k], min_size * msa); coil, pfcoil, r0, z0, width0, height0)
        end
    end
end

#= =================================================================================== =#
#  Visualization of IMAS.pf_active.coil, pf_active__supply, pf_passive__loop as tables  #
#= =================================================================================== =#
function Base.show(io::IO, mime::MIME"text/plain", coils::IMAS.IDSvector{<:Union{IMAS.pf_active__coil,IMAS.pf_active__supply,IMAS.pf_passive__loop}})
    old_lines = get(ENV, "LINES", missing)
    old_columns = get(ENV, "COLUMNS", missing)
    df = DataFrames.DataFrame(coils)
    try
        ENV["LINES"] = 1000
        ENV["COLUMNS"] = 1000
        return show(io, mime, df)
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

function DataFrames.DataFrame(loops::IMAS.IDSvector{<:IMAS.pf_passive__loop{T}}) where {T<:Real}

    df = DataFrames.DataFrame(;
        name=String[],
        n_elements=Int[],
        n_total_turns=T[],
        resistivity=T[]
    )

    for loop in loops
        turns = sum(getproperty(element, :turns_with_sign, 1.0) for element in loop.element)
        resistivity = getproperty(loop, :resistivity, NaN)
        push!(df, [loop.name, length(loop.element), turns, resistivity])
    end

    return df
end

function DataFrames.DataFrame(coils::IMAS.IDSvector{<:IMAS.pf_active__coil{T}}) where {T<:Real}

    df = DataFrames.DataFrame(;
        name=String[],
        var"function"=Vector{Symbol}[],
        n_elements=Int[],
        n_total_turns=T[],
        resitance=T[]
    )

    for coil in coils
        func = [IMAS.index_2_name(coil.function)[f.index] for f in coil.function]
        turns = sum(getproperty(element, :turns_with_sign, 1.0) for element in coil.element)
        resitance = getproperty(coil, :resistance, NaN)
        push!(df, [coil.name, func, length(coil.element), sum(turns), resitance])
    end

    return df
end

function DataFrames.DataFrame(supplies::IMAS.IDSvector{<:IMAS.pf_active__supply})

    df = DataFrames.DataFrame(;
        name=String[],
        voltage_limits=Tuple{Float64,Float64}[],
        current_limits=Tuple{Float64,Float64}[]
    )

    for supply in supplies
        push!(df, [supply.name, (supply.voltage_limit_min, supply.voltage_limit_max), (supply.current_limit_min, supply.current_limit_max)])
    end

    return df
end

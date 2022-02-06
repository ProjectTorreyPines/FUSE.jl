using Equilibrium
import PolygonOps
using StaticArrays
using Interpolations
using Optim
using AD_GS
using LinearAlgebra
using Statistics
using Plots
import Contour

const coils_turns_spacing = 0.05

#= ================== =#
#  init pf_active IDS  #
#= ================== =#
function finite_size_OH_coils(z, coils_cleareance)
    ez = diff(z) / 2.0 .+ z[1:end-1]
    ez = vcat((ez[1] - ez[2]) + ez[1], ez, (ez[end] - ez[end-1]) + ez[end])
    ez = (ez .- minimum(ez)) ./ (maximum(ez) - minimum(ez)) * (maximum(z) - minimum(z)) .+ minimum(z)
    ez_centers = diff(ez) / 2.0 .+ ez[1:end-1]
    ez_centers = [abs(z)<1E-6 ? 0 : z for z in ez_centers] # correct small deviations near zero
    ez_heights = diff(ez) .- coils_cleareance
    return ez_centers, ez_heights
end

"""
    init(pf_active::IMAS.pf_active, bd::IMAS.build, n_coils::Vector)

Use build layers outline to initialize PF coils distribution

Attributes
 * n_coils: number of pf coils per coil-placement rail (the first one is the OH, and then one for each vacuum region)
 * pf_coils_size:  Size of the (square) coils (per PF rail)
 * coils_cleareance: Clereance that coils have from other structures (per rail)
 * coils_elements_area: Cross-sectional area taken up by individual filaments in a coil (per rail)
"""
function init(pf_active::IMAS.pf_active,
              bd::IMAS.build,
              n_coils::Vector{TI} where TI <: Int;
              pf_coils_size::Union{Nothing, Real, Vector{TR}} where TR <: Real=nothing,
              coils_cleareance::Union{Nothing, Real, Vector{TR}} where TR <: Real=nothing,
              coils_elements_area::Union{Nothing, Real, Vector{TR}} where TR <: Real=nothing)

    OH_layer = IMAS.get_build(bd, type=1)

    empty!(pf_active)
    resize!(bd.pf_coils_rail, length(n_coils))

    # make sure coils_cleareance is an array the lenght of the rails
    if coils_elements_area === nothing
        coils_elements_area = 0.0025
    end
    if isa(coils_elements_area, Number)
        coils_elements_area = [coils_elements_area for k in 1:length(n_coils)]
    end

    # make sure coils_cleareance is an array the lenght of the rails
    if coils_cleareance === nothing
        coils_cleareance = (maximum(OH_layer.outline.r) - minimum(OH_layer.outline.r)) / 2.0
    end
    if isa(coils_cleareance, Number)
        coils_cleareance = [coils_cleareance for k in 1:length(n_coils)]
    end

    # OH coils are distributed on a rail within the OH region
    r_ohcoils = ones(n_coils[1]) .* (sum(extrema(OH_layer.outline.r)) / 2.)
    w = maximum(OH_layer.outline.r) - minimum(OH_layer.outline.r)
    z_ohcoils = collect(range(minimum(OH_layer.outline.z), maximum(OH_layer.outline.z), length=n_coils[1]))
    z_ohcoils, h_ohcoils = finite_size_OH_coils(z_ohcoils, coils_cleareance[1])
    bd.pf_coils_rail[1].name = "OH"
    bd.pf_coils_rail[1].coils_number = n_coils[1]
    bd.pf_coils_rail[1].coils_elements_area = coils_elements_area[1]
    bd.pf_coils_rail[1].coils_cleareance = coils_cleareance[1]
    bd.pf_coils_rail[1].outline.r = r_ohcoils
    bd.pf_coils_rail[1].outline.z = z_ohcoils
    bd.pf_coils_rail[1].outline.distance = range(-1, 1, length=n_coils[1])
    for (r, z, h) in zip(r_ohcoils, z_ohcoils, h_ohcoils)
        k = length(pf_active.coil) + 1
        resize!(pf_active.coil, k)
        resize!(pf_active.coil[k].element, 1)
        pf_active.coil[k].identifier = "optim"
        pf_active.coil[k].name = "OH"
        pf_active.coil[k].element[1].geometry.rectangle.r = r
        pf_active.coil[k].element[1].geometry.rectangle.z = z
        pf_active.coil[k].element[1].geometry.rectangle.width = w
        pf_active.coil[k].element[1].geometry.rectangle.height = h
        set_turns_from_spacing!(pf_active.coil[k], coils_turns_spacing, +1)
        @ddtime pf_active.coil[k].current.data = 0.0
    end

    # make sure coils_cleareance is an array the lenght of the PF rails
    if pf_coils_size === nothing
        pf_coils_size = sqrt(w * sum(h_ohcoils) / length(h_ohcoils))
    end
    if isa(pf_coils_size, Number)
        pf_coils_size = reverse([pf_coils_size/s for s in 2.0.^(0:length(n_coils)-2)])
    end

    # Now add actual PF coils to regions of vacuum
    krail = 1
    resolution = 257
    rmask, zmask, mask = IMAS.structures_mask(bd, resolution=resolution)
    for (k, layer) in enumerate(bd.layer)
        if (layer.hfs == 1 || k == length(bd.layer)) && ! ismissing(layer.outline, :r)
            if ! ismissing(layer, :material) && layer.material == "vacuum"

                krail += 1
                nc = n_coils[krail]

                # add rail info to build IDS
                bd.pf_coils_rail[krail].name = replace(replace(layer.name, "hfs " => ""), "lfs " => "")
                bd.pf_coils_rail[krail].coils_number = nc
                bd.pf_coils_rail[krail].coils_elements_area = coils_elements_area[krail]
                bd.pf_coils_rail[krail].coils_cleareance = coils_cleareance[krail]

                # pick layers with outline information
                if layer.hfs == 1
                    outer_layer = IMAS.get_build(bd, identifier=bd.layer[k].identifier, hfs=1)
                    inner_layer = IMAS.get_build(bd, identifier=bd.layer[k + 1].identifier, hfs=[1,0])
                else
                    inner_layer = IMAS.get_build(bd, identifier=bd.layer[k - 1].identifier, hfs=1)
                    outer_layer = IMAS.get_build(bd, identifier=bd.layer[k].identifier, hfs=[1,0])
                end

                # generate rail between the two layers where coils will be placed and will be able to slide during the `optimization` phase
                coil_size = pf_coils_size[krail - 1]
                poly = LibGEOS.buffer(xy_polygon(inner_layer.outline.r, inner_layer.outline.z), coil_size / sqrt(2) + coils_cleareance[krail])
                rail_r = [v[1] for v in LibGEOS.coordinates(poly)[1]]
                rail_z = [v[2] for v in LibGEOS.coordinates(poly)[1]]

                # mark what regions on that rail do not intersect solid structures and can hold coils
                iclearance = Int(floor(coil_size/(rmask[2] - rmask[1])/2))
                valid_k = []
                for (k, (r, z)) in enumerate(zip(rail_r, rail_z))
                    ir = argmin(abs.(rmask .- r))
                    iz = argmin(abs.(zmask .- z))
                    if (ir - iclearance) < 1 || (ir + iclearance) > length(rmask) || (iz - iclearance) < 1 || (iz + iclearance) > length(zmask)
                        continue
                    end
                    if all(mask[(-iclearance:iclearance) .+ ir,(-iclearance:iclearance) .+ iz] .== 0)
                        push!(valid_k, k)
                    end
                end
                if length(valid_k) == 0
                    bd.pf_coils_rail[krail].outline.r = Real[]
                    bd.pf_coils_rail[krail].outline.z = Real[]
                    bd.pf_coils_rail[krail].outline.distance = Real[]
                    error("Coils on PF rail #$(krail-1) are too big to fit.")
                    continue
                end
                istart = argmax(diff(valid_k)) + 1
                if istart < (valid_k[1] + (length(rail_r) - valid_k[end]))
                    istart = 0
                end
                valid_r = fill(NaN, size(rail_r)...)
                valid_z = fill(NaN, size(rail_z)...)
                valid_r[valid_k] = rail_r[valid_k]
                valid_z[valid_k] = rail_z[valid_k]
                valid_r = vcat(valid_r[istart + 1:end], valid_r[1:istart])
                valid_z = vcat(valid_z[istart + 1:end], valid_z[1:istart])

                # evaluate distance along rail
                d_distance = sqrt.(diff(vcat(valid_r, valid_r[1])).^2.0 .+ diff(vcat(valid_z, valid_z[1])).^2.0)
                d_distance[isnan.(d_distance)] .= 0.0
                distance = cumsum(d_distance)
                valid_z = valid_z[d_distance .!= 0]
                valid_r = valid_r[d_distance .!= 0]
                distance = distance[d_distance .!= 0]
                distance = (distance .- distance[1])
                distance = (distance ./ distance[end]).*2.0.-1.0

                # add rail info to build IDS
                bd.pf_coils_rail[krail].outline.r = valid_r
                bd.pf_coils_rail[krail].outline.z = valid_z
                bd.pf_coils_rail[krail].outline.distance = distance

                if nc == 0
                    continue
                end

                # uniformely distribute coils
                coils_distance = range(-(1-1/nc),1-1/nc,length=nc)
                r_coils = IMAS.interp(distance, valid_r)(coils_distance)
                z_coils = IMAS.interp(distance, valid_z)(coils_distance)
                z_coils = [abs(z)<1E-6 ? 0 : z for z in z_coils]

                # populate IMAS data structure
                for (r, z) in zip(r_coils, z_coils)
                    k = length(pf_active.coil) + 1
                    resize!(pf_active.coil, k)
                    resize!(pf_active.coil[k].element, 1)
                    pf_active.coil[k].identifier = "optim"
                    pf_active.coil[k].name = "PF"
                    pf_active.coil[k].element[1].geometry.rectangle.r = r
                    pf_active.coil[k].element[1].geometry.rectangle.z = z
                    pf_active.coil[k].element[1].geometry.rectangle.width = coil_size
                    pf_active.coil[k].element[1].geometry.rectangle.height = coil_size
                    set_turns_from_spacing!(pf_active.coil[k], coils_turns_spacing, +1)
                    @ddtime pf_active.coil[k].current.data = 0.0
                end
            end
        end
    end

    return pf_active
end

#= =============== =#
#  PFcoilsOptActor  #
#= =============== =#
@Base.kwdef mutable struct PFcoilsOptTrace
    params::Vector{Vector{Real}} = Vector{Real}[]
    cost_ψ::Vector{Real} = Real[]
    cost_currents::Vector{Real} = Real[]
    cost_total::Vector{Real} = Real[]
end

mutable struct PFcoilsOptActor <: AbstractActor
    eq_in::IMAS.equilibrium
    eq_out::IMAS.equilibrium
    time::Real
    pf_active::IMAS.pf_active
    bd::IMAS.build
    symmetric::Bool
    λ_regularize::Real
    trace::PFcoilsOptTrace
    coil_model::Symbol
end

function PFcoilsOptActor(eq_in::IMAS.equilibrium,
                         bd::IMAS.build,
                         pf::IMAS.pf_active,
                         n_coils::Vector;
                         λ_regularize=1E-13,
                         coil_model=:simple)
    # initialize coils location
    init(pf, bd, n_coils)

    # basic constructors
    eq_out = deepcopy(eq_in)
    symmetric = false
    time_index = 1
    time = eq_in.time[time_index]

    # constructor
    pfactor = PFcoilsOptActor(eq_in, eq_out, time, pf, bd, symmetric, λ_regularize, PFcoilsOptTrace(), coil_model)

    return pfactor
end

# Dispatching AD_GS on IMAS.pf_active__coil
mutable struct GS_IMAS_pf_active__coil <: AD_GS.AbstractCoil
    pf_active__coil::IMAS.pf_active__coil
    r::Real
    z::Real
    width::Real
    height::Real
    turns_with_sign::Real
    spacing::Real
    time_current::Vector{T} where T <: Real
    time::Vector{T} where T <: Real
    time_index::Int
    coil_model::Symbol
end

function GS_IMAS_pf_active__coil(pf_active__coil, coil_model)
    return GS_IMAS_pf_active__coil(pf_active__coil,
                                   pf_active__coil.element[1].geometry.rectangle.r,
                                   pf_active__coil.element[1].geometry.rectangle.z,
                                   pf_active__coil.element[1].geometry.rectangle.width,
                                   pf_active__coil.element[1].geometry.rectangle.height,
                                   pf_active__coil.element[1].turns_with_sign,
                                   get_spacing_from_turns(pf_active__coil),
                                   pf_active__coil.current.data,
                                   pf_active__coil.current.time,
                                   1,
                                   coil_model)
end

function Base.getproperty(coil::GS_IMAS_pf_active__coil, field::Symbol)
    if field == :current
        return getfield(coil,:time_current)[coil.time_index]
    else
        return getfield(coil, field)
    end
end

function Base.setproperty!(coil::GS_IMAS_pf_active__coil, field::Symbol, value)
    if field == :current
        getfield(coil,:time_current)[coil.time_index] = value
    else
        setfield!(coil, field, value)
    end
    if field in [:width, :height, :spacing]
        s = sign(getfield(coil, :turns_with_sign))
        turns = Int(ceil(coil.width .* coil.height ./ coil.spacing.^2))
        setfield!(coil, :turns_with_sign, s * turns)
    end
end

function transfer_info_GS_coil_to_IMAS(coil::GS_IMAS_pf_active__coil)
    pf_active__coil = coil.pf_active__coil
    pf_active__coil.element[1].geometry.rectangle.r = coil.r
    pf_active__coil.element[1].geometry.rectangle.z = coil.z
    pf_active__coil.element[1].geometry.rectangle.width = coil.width
    pf_active__coil.element[1].geometry.rectangle.height = coil.height
    pf_active__coil.element[1].turns_with_sign = coil.turns_with_sign
    pf_active__coil.current.time = coil.time
    pf_active__coil.current.data = coil.time_current
end

function set_turns_from_spacing!(coil::GS_IMAS_pf_active__coil)
    pf_active__coil = getfield(coil,:pf_active__coil)
    return set_turns_from_spacing!(pf_active__coil, coil.spacing)
end

function set_turns_from_spacing!(pf_active__coil::IMAS.pf_active__coil, spacing::Real)
    s = sign(pf_active__coil.element[1].turns_with_sign)
    set_turns_from_spacing!(pf_active__coil, spacing, s)
end

function set_turns_from_spacing!(pf_active__coil::IMAS.pf_active__coil, spacing::Real, s::Int)
    area = (pf_active__coil.element[1].geometry.rectangle.width * pf_active__coil.element[1].geometry.rectangle.height)
    pf_active__coil.element[1].turns_with_sign = s * Int(ceil(area / spacing^2))
end

function get_spacing_from_turns(coil::GS_IMAS_pf_active__coil)
    pf_active__coil = getfield(coil,:pf_active__coil)
    return get_spacing_from_turns(pf_active__coil)
end

function get_spacing_from_turns(pf_active__coil::IMAS.pf_active__coil)
    return sqrt((pf_active__coil.element[1].geometry.rectangle.width * pf_active__coil.element[1].geometry.rectangle.height) / abs(pf_active__coil.element[1].turns_with_sign))
end

"""
    AD_GS.Green(coil::GS_IMAS_pf_active__coil, R::Real, Z::Real)

Calculates coil green function at given R and Z coordinate
"""
function AD_GS.Green(coil::GS_IMAS_pf_active__coil, R::Real, Z::Real)
    # fastest
    if coil.coil_model == :point
        return AD_GS.Green(coil.r, coil.z, R, Z, coil.turns_with_sign)

    # medium
    elseif coil.coil_model in [:corners, :simple]
        if coil.pf_active__coil.name == "OH"
            n = 3
            z_filaments = range(coil.z - (coil.height - coil.width/2.0) / 2.0, coil.z + (coil.height - coil.width/2.0) / 2.0, length=n)
            green = []
            for z in z_filaments
                push!(green, AD_GS.Green(coil.r, z, R, Z, coil.turns_with_sign / n))
            end
            return sum(green)
        elseif coil.coil_model == :corners
            return AD_GS.Green(AD_GS.ParallelogramCoil(coil.r, coil.z, coil.width/2.0, coil.height/2.0, 0.0, 90.0, nothing), R, Z, coil.turns_with_sign/4)
        elseif coil.coil_model == :simple
            return AD_GS.Green(coil.r, coil.z, R, Z, coil.turns_with_sign)
        end

    # slow
    elseif coil.coil_model == :realistic
        return AD_GS.Green(AD_GS.ParallelogramCoil(coil.r, coil.z, coil.width, coil.height, 0.0, 90.0, coil.spacing), R, Z)

    else
        error("GS_IMAS_pf_active__coil coil.coil_model can only be (in order of accuracy) :realistic, :corners, :simple, and :point")
    end
end

# step
function pack_rail(bd::IMAS.build, λ_regularize::Float64, symmetric::Bool)::Vector{Float64}
    distances = []
    for rail in bd.pf_coils_rail
        if rail.name !== "OH"
            # not symmetric
            if ! symmetric
                coil_distances = collect(range(-1.0, 1.0, length=rail.coils_number + 2))[2:end - 1]
            # even symmetric
            elseif mod(rail.coils_number, 2) == 0
                coil_distances = collect(range(-1.0, 1.0, length=rail.coils_number + 2))[2 + Int(rail.coils_number // 2):end - 1]
            # odd symmetric
            else
                coil_distances = collect(range(-1.0, 1.0, length=rail.coils_number + 2))[2 + Int((rail.coils_number - 1) // 2) + 1:end - 1]
            end
            append!(distances, coil_distances)
        end
    end
    oh_height_off=[]
    for rail in bd.pf_coils_rail
        if rail.name == "OH"
            push!(oh_height_off, 1.0)
            if ! symmetric
                push!(oh_height_off, 0.0)
            end
        end
    end
    packed = vcat(distances, oh_height_off, log10(λ_regularize))

    return packed
end

function unpack_rail!(packed::Vector, optim_coils::Vector, symmetric::Bool, bd::IMAS.build)
    λ_regularize = packed[end]
    if symmetric
        n_oh_params=1
    else
        n_oh_params=2
    end
    if any(rail.name == "OH" for rail in bd.pf_coils_rail)
        oh_height_off = packed[end - n_oh_params:end - 1]
        distances = packed[1:end - n_oh_params]
    else
        oh_height_off = []
        distances = packed[1:end - 1]
    end

    if length(optim_coils) != 0 # optim_coils have zero length in case of the `static` optimization
        kcoil = 0
        koptim = 0
        koh = 0
        for rail in bd.pf_coils_rail
            if rail.name == "OH"
                for k in 1:rail.coils_number
                    koptim += 1
                    koh += 1

                    # mirror OH size when it reaches maximum extent of the rail
                    while (oh_height_off[1] < -1) || (oh_height_off[1] > 1)
                        if oh_height_off[1] < -1
                            oh_height_off[1] = -2.0 .- oh_height_off[1]
                        else
                            oh_height_off[1] = 2.0 .- oh_height_off[1]
                        end
                    end
                    Δrail = maximum(rail.outline.z)-minimum(rail.outline.z)
                    rail_offset = (maximum(rail.outline.z)+minimum(rail.outline.z))/2.0
                    optim_coils[koptim].z = range(-oh_height_off[1]/2.0,oh_height_off[1]/2.0,length=rail.coils_number)[koh]*Δrail + rail_offset
                    if ! symmetric
                        optim_coils[koptim].z += oh_height_off[2] * (1-oh_height_off[1]) * Δrail
                    end
                    optim_coils[koptim].height = (oh_height_off[1] * Δrail - rail.coils_cleareance * (rail.coils_number - 1)) / rail.coils_number
                end
            else
                r_interp = IMAS.interp(rail.outline.distance, rail.outline.r, extrapolation_bc=:flat)
                z_interp = IMAS.interp(rail.outline.distance, rail.outline.z, extrapolation_bc=:flat)
                # not symmetric
                if ! symmetric
                    dkcoil = rail.coils_number
                    coil_distances = distances[kcoil + 1:kcoil + dkcoil]
                # even symmetric
                elseif mod(rail.coils_number, 2) == 0
                    dkcoil = Int(rail.coils_number // 2)
                    coil_distances = distances[kcoil + 1:kcoil + dkcoil]
                    coil_distances = vcat(- reverse(coil_distances), coil_distances)
                # odd symmetric
                else
                    dkcoil = Int((rail.coils_number - 1) // 2)
                    coil_distances = distances[kcoil + 1:kcoil + dkcoil]
                    coil_distances = vcat(- reverse(coil_distances), 0.0, coil_distances)
                end
                kcoil += dkcoil

                # mirror coil position when they reach the end of the rail
                while any(coil_distances .< -1) || any(coil_distances .> 1)
                    coil_distances[coil_distances .< -1] = -2.0 .- coil_distances[coil_distances .< -1]
                    coil_distances[coil_distances .> 1] = 2.0 .- coil_distances[coil_distances .> 1]
                end

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

function optimize_coils_rail(eq::IMAS.equilibrium;pinned_coils::Vector, optim_coils::Vector, fixed_coils::Vector, symmetric::Bool, λ_regularize::Real, λ_ψ::Real, λ_null::Real, λ_currents::Real, λ_strike::Real, bd::IMAS.build, maxiter::Int, verbose::Bool)
    fixed_eqs = []
    weights = []
    for time_index in 1:length(eq.time_slice)
        eqt = eq.time_slice[time_index]
        # field nulls
        if ismissing(eqt.global_quantities, :ip)
            # find ψp
            Bp_fac, ψp, Rp, Zp = AD_GS.field_null_on_boundary(eqt.global_quantities.psi_boundary,
                                                              eqt.boundary.outline.r,
                                                              eqt.boundary.outline.z,
                                                              fixed_coils)
            push!(fixed_eqs, (Bp_fac, ψp, Rp, Zp))
            push!(weights, nothing)
        # solutions with plasma
        else
            fixed_eq = IMAS2Equilibrium(eqt)
            # private flux regions
            private = IMAS.flux_surface(eqt,eqt.profiles_1d.psi[end],false)
            vessel = IMAS.get_build(bd, type=-1, hfs=0)
            Rx = []
            Zx = []
            for (pr, pz) in private
               pvx , pvy = IMAS.intersection(vessel.outline.r, vessel.outline.z, pr, pz; as_list_of_points=false)
               append!(Rx, pvx)
               append!(Zx, pvy)
            end
            # find ψp
            Bp_fac, ψp, Rp, Zp = AD_GS.ψp_on_fixed_eq_boundary(fixed_eq, fixed_coils; Rx, Zx)
            push!(fixed_eqs, (Bp_fac, ψp, Rp, Zp))
            # give each strike point the same weight as the lcfs
            weight = Rp .* 0.0 .+ 1.0
            weight[end - length(Rx) + 1 : end] .= length(Rp) / (1+length(Rx)) * λ_strike
            push!(weights, weight)
        end
    end

    packed = pack_rail(bd, λ_regularize, symmetric)
    trace = PFcoilsOptTrace()

    packed_tmp = []
    function placement_cost(packed; do_trace=false)
        try
            push!(packed_tmp, packed)
            λ_regularize = unpack_rail!(packed, optim_coils, symmetric, bd)
            coils = vcat(pinned_coils, optim_coils)
            all_cost_ψ = []
            all_cost_currents = []
            for (time_index, (fixed_eq, weight)) in enumerate(zip(fixed_eqs,weights))
                for coil in vcat(pinned_coils, optim_coils, fixed_coils)
                    coil.time_index = time_index
                end
                currents, cost_ψ0 = AD_GS.currents_to_match_ψp(fixed_eq..., coils, weights=weight, λ_regularize=λ_regularize, return_cost=true)
                if ismissing(eq.time_slice[time_index].global_quantities, :ip)
                    push!(all_cost_ψ, cost_ψ0 / λ_null)
                else
                    push!(all_cost_ψ, cost_ψ0 / λ_ψ)
                end
                push!(all_cost_currents, norm((exp.(currents/λ_currents).-1.0)/(exp(1)-1)) / length(currents))
            end
            cost_ψ = norm(all_cost_ψ) / length(all_cost_ψ)
            cost_currents = norm(all_cost_currents) / length(all_cost_currents)
            cost = sqrt(cost_ψ^2 + cost_currents^2)
            if do_trace
                push!(trace.params, packed)
                push!(trace.cost_ψ, cost_ψ)
                push!(trace.cost_currents, cost_currents)
                push!(trace.cost_total, cost)
            end
            return cost
        catch e
            println(e)
            rethrow
        end

    end

    function clb(x)
        placement_cost(packed_tmp[end]; do_trace=true)
        false
    end
    
    if maxiter == 0
        placement_cost(packed)
        λ_regularize = unpack_rail!(packed, optim_coils, symmetric, bd)
    else
        # use NelderMead() ; other optimizer that works is Newton(), others have trouble
        res = Optim.optimize(placement_cost, packed, Optim.NelderMead(), Optim.Options(time_limit=60 * 2, iterations=maxiter, callback=clb); autodiff=:forward)
        if verbose println(res) end
        packed = Optim.minimizer(res)
        λ_regularize = unpack_rail!(packed, optim_coils, symmetric, bd)
    end

    return λ_regularize, trace
end


"""
    fixed_pinned_optim_coils(pfactor, optimization_scheme)

Returns tuple of GS_IMAS_pf_active__coil coils organized by their function:
- fixed: fixed position and current
- pinned: coisl with fixed position but current is optimized
- optim: coils that have theri position and current optimized
"""
function fixed_pinned_optim_coils(pfactor, optimization_scheme)
    fixed_coils = GS_IMAS_pf_active__coil[]
    pinned_coils = GS_IMAS_pf_active__coil[]
    optim_coils = GS_IMAS_pf_active__coil[]
    for coil in pfactor.pf_active.coil
        if coil.identifier == "pinned"
            push!(pinned_coils, GS_IMAS_pf_active__coil(coil, pfactor.coil_model))
        elseif (coil.identifier == "optim") && (optimization_scheme == :static)
            push!(pinned_coils, GS_IMAS_pf_active__coil(coil, pfactor.coil_model))
        elseif coil.identifier == "optim"
            push!(optim_coils, GS_IMAS_pf_active__coil(coil, pfactor.coil_model))
        elseif coil.identifier == "fixed"
            push!(fixed_coils, GS_IMAS_pf_active__coil(coil, pfactor.coil_model))
        else
            error("Accepted type of coil.identifier are only \"optim\", \"pinned\", or \"fixed\"")
        end
    end
    return fixed_coils, pinned_coils, optim_coils
end

"""
    step(pfactor::PFcoilsOptActor;
        symmetric=pfactor.symmetric,
        λ_regularize=pfactor.λ_regularize,
        λ_ψ=1E-2,
        λ_null=1,
        λ_currents=1E5,
        λ_strike=1,
        maxiter=10000,
        optimization_scheme=:rail,
        verbose=false)

Optimize coil currents and positions to produce sets of equilibria while minimizing coil currents
"""
function step(pfactor::PFcoilsOptActor;
              symmetric=pfactor.symmetric,
              λ_regularize=pfactor.λ_regularize,
              λ_ψ=1E-2,
              λ_null=1,
              λ_currents=1E5,
              λ_strike=1,
              maxiter=10000,
              optimization_scheme=:rail,
              verbose=false)

    fixed_coils, pinned_coils, optim_coils = fixed_pinned_optim_coils(pfactor, optimization_scheme)
    coils = vcat(pinned_coils, optim_coils, fixed_coils)
    for coil in coils
        coil.time_current = pfactor.eq_in.time .* 0.0
        coil.time = pfactor.eq_in.time
    end

    bd = pfactor.bd
    # run rail type optimizer
    if optimization_scheme in [:rail, :static]
        (λ_regularize, trace) = optimize_coils_rail(pfactor.eq_in; pinned_coils, optim_coils, fixed_coils, symmetric, λ_regularize, λ_ψ, λ_null, λ_currents, λ_strike, bd, maxiter, verbose)
    else
        error("Supported PFcoilsOptActor optimization_scheme are `:static` or `:rail`")
    end
    pfactor.λ_regularize = λ_regularize
    pfactor.trace = trace

    # transfer the results to IMAS.pf_active
    for coil in coils
        transfer_info_GS_coil_to_IMAS(coil)
    end

    return pfactor
end


"""
    finalize(pfactor::PFcoilsOptActor; scale_eq_domain_size = 1.0)

Update pfactor.eq_out 2D equilibrium PSI based on coils positions and currents
"""
function finalize(pfactor::PFcoilsOptActor; scale_eq_domain_size = 1.0)
    coils = GS_IMAS_pf_active__coil[]
    for coil in pfactor.pf_active.coil
        push!(coils, GS_IMAS_pf_active__coil(coil, pfactor.coil_model))
    end

    # update equilibrium
    for time_index in 1:length(pfactor.eq_in.time_slice)
        if ismissing(pfactor.eq_in.time_slice[time_index].global_quantities, :ip)
            continue
        end
        for coil in coils
            coil.time_index = time_index
        end

        # convert equilibrium to Equilibrium.jl format, since this is what AD_GS uses
        EQfixed = IMAS2Equilibrium(pfactor.eq_in.time_slice[time_index])

        # # update ψ map
        R = range(EQfixed.r[1] / scale_eq_domain_size, EQfixed.r[end] * scale_eq_domain_size, length = length(EQfixed.r))
        Z = range(EQfixed.z[1] * scale_eq_domain_size, EQfixed.z[end] * scale_eq_domain_size, length = length(EQfixed.z))
        ψ_f2f = AD_GS.fixed2free(EQfixed, coils, R, Z)
        pfactor.eq_out.time_slice[time_index].profiles_2d[1].grid.dim1 = R
        pfactor.eq_out.time_slice[time_index].profiles_2d[1].grid.dim2 = Z
        pfactor.eq_out.time_slice[time_index].profiles_2d[1].psi = transpose(ψ_f2f)
        # IMAS.flux_surfaces(pfactor.eq_out.time_slice[time_index]) #### PROBLEM
    end
end

# plotting
"""
    plot_pfcoilsactor_cx(pfactor::PFcoilsOptActor; time_index=1, equilibrium=true, rail=true)

Plot PFcoilsOptActor optimization cross-section
"""
@recipe function plot_pfcoilsactor_cx(pfactor::PFcoilsOptActor; time_index=1, equilibrium=true, build=true, coils_flux=false, rail=false, plot_r_buffer=1.6)

    # if there is no equilibrium then treat this as a field_null plot
    field_null = false
    if length(pfactor.eq_out.time_slice[time_index].profiles_2d)==0 || ismissing(pfactor.eq_out.time_slice[time_index].profiles_2d[1], :psi)
        coils_flux = true
        field_null = true
    end

    # when plotting coils_flux the build is not visible anyways
    if coils_flux
        build = false
    end

    # setup plotting area
    xlim = [0.0, maximum(pfactor.bd.layer[end].outline.r)]
    ylim = [minimum(pfactor.bd.layer[end].outline.z), maximum(pfactor.bd.layer[end].outline.z)]
    xlim --> xlim * plot_r_buffer
    ylim --> ylim
    aspect_ratio --> :equal

    # plot build
    if build
        @series begin
            exclude_layers --> [:oh]
            pfactor.bd
        end
    end

    # plot coils_flux
    if coils_flux
        resolution = 129
        R = range(xlim[1], xlim[2], length=resolution)
        Z = range(ylim[1], ylim[2], length=Int(ceil(resolution * (ylim[2]-ylim[1]) / (xlim[2]- xlim[1]))))

        coils = [GS_IMAS_pf_active__coil(coil, pfactor.coil_model) for coil in pfactor.pf_active.coil]
        for coil in coils
            coil.time_index=time_index
        end

        # ψ coil currents
        ψbound = pfactor.eq_out.time_slice[time_index].global_quantities.psi_boundary
        ψ = AD_GS.coils_flux(2*pi, coils, R, Z)

        ψmin = minimum(x->isnan(x) ? Inf : x, ψ)
        ψmax = maximum(x->isnan(x) ? -Inf : x, ψ)
        ψabsmax = maximum(x->isnan(x) ? -Inf : x, abs.(ψ))
        
        if field_null
            clims = (-ψabsmax / 10 + ψbound, ψabsmax / 10 + ψbound)
        else
            clims = (ψmin, ψmax)
        end

        @series begin
            seriestype --> :contourf
            c --> :diverging
            colorbar_entry --> false
            levels --> range(clims[1],clims[2],length=21)
            linewidth --> 0.0
            R, Z, transpose(ψ)
        end

        if field_null
            @series begin
                seriestype --> :contour
                colorbar_entry --> false
                levels --> [ψbound]
                linecolor --> :black
                R, Z, transpose(ψ)
            end
        end

        @series begin
            outlines --> true
            exclude_layers --> [:oh]
            pfactor.bd
        end
    end

    # plot equilibrium
    if equilibrium
        if field_null
            @series begin
                label --> "Field null region"
                seriescolor --> :red
                pfactor.eq_out.time_slice[time_index]
            end
        else
            @series begin
                label --> "Final"
                seriescolor --> :red
                pfactor.eq_out.time_slice[time_index]
            end
            @series begin
                label --> "Target"
                seriescolor --> :blue
                lcfs --> true
                linestyle --> :dash
                pfactor.eq_in.time_slice[time_index]
            end
        end
    end

    # plot pf_active coils
    @series begin
        time_index --> time_index
        pfactor.pf_active
    end

    # plot optimization rails
    if rail
        for (krail, rail) in enumerate(pfactor.bd.pf_coils_rail)
            if ! ismissing(rail.outline,:r)
                @series begin
                    label --> (build ? "Coil opt. rail" : "")
                    primary --> krail == 1 ? true : false
                    color --> :gray
                    linestyle --> :dash
                    rail.outline.r, rail.outline.z
                end
            end
        end
    end

end

"""
    plot_pfcoilsactor_trace(trace::PFcoilsOptTrace, what::Symbol=:cost; start_at::Int=1)

Plot PFcoilsOptActor optimization trace

Attributes:
- what::Symbol=:cost or :currents or individual fields of the PFcoilsOptTrace structure
- start_at=::Int=1 index of the first element of the trace to start plotting
"""
@recipe function plot_pfcoilsactor_trace(trace::PFcoilsOptTrace, what::Symbol=:cost; start_at=1)
    start_at = minimum([start_at, length(trace.cost_total)])
    x = start_at:length(trace.cost_total)
    legend --> :bottomleft
    if what == :cost
        if sum(trace.cost_ψ[start_at:end]) > 0.0
            @series begin
                label --> "ψ"
                yscale --> :log10
                x, trace.cost_ψ[start_at:end]
            end
        end
        if sum(trace.cost_currents[start_at:end]) > 0.0
            @series begin
                label --> "currents"
                yscale --> :log10
                x, trace.cost_currents[start_at:end]
            end
        end
        @series begin
            label --> "total"
            yscale --> :log10
            linestyle --> :dash
            color --> :black
            # ylim --> [minimum(trace.cost_total[start_at:end]) / 10,maximum(trace.cost_total[start_at:end])]
            x, trace.cost_total[start_at:end]
        end

    elseif what == :currents
        @series begin
            label --> "Starting"
            getfield(trace, what)[start_at:end][1,:]
        end
        @series begin
            label --> "Final"
            getfield(trace, what)[start_at:end][end,:]
        end

    else
        @series begin
            if occursin("cost_", String(what))
                yscale --> :log10
            end
            label --> String(what)
            x, getfield(trace, what)[start_at:end]
        end
    end
end

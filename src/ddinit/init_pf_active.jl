#= ================== =#
#  init pf_active IDS  #
#= ================== =#
"""
    init_pf_active!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())

Initialize `dd.pf_active` starting from `ini` and `act` parameters
"""
function init_pf_active!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())
    TimerOutputs.reset_timer!("init_pf_active")
    TimerOutputs.@timeit timer "init_pf_active" begin
        init_from = ini.general.init_from
        if init_from == :ods
            if length(dd1.pf_active.coil) > 0
                dd.pf_active = deepcopy(dd1.pf_active)
                IMAS.set_coils_function(dd.pf_active.coil, ini.equilibrium.R0)
            else
                init_from = :scalars
            end
        end

        if init_from == :scalars
            init_pf_active!(dd.pf_active, dd.build, dd.equilibrium.time_slice[])
        end

        IMAS.coil_technology(dd.build.tf.technology, ini.tf.technology, :tf)
        IMAS.coil_technology(dd.build.oh.technology, ini.oh.technology, :oh)
        IMAS.coil_technology(dd.build.pf_active.technology, ini.pf_active.technology, :pf_active)

        return dd
    end
end

function clip_rails(
    rail_r::AbstractVector{T1},
    rail_z::AbstractVector{T1},
    pr::AbstractVector{T2},
    pz::AbstractVector{T2},
    RA::T3,
    ZA::T3,
    side::IMAS.BuildLayerSide
) where {T1<:Real,T2<:Real,T3<:Real}
    rail_r = convert(Vector{Float64}, rail_r)
    rail_z = convert(Vector{Float64}, rail_z)
    pr = convert(Vector{Float64}, pr)
    pz = convert(Vector{Float64}, pz)
    RA = convert(Float64, RA)
    ZA = convert(Float64, ZA)
    return clip_rails(rail_r, rail_z, pr, pz, RA, ZA, side)
end

"""
    clip_rails(rail_r::Vector{T}, rail_z::Vector{T}, pr::Vector{T}, pz::Vector{T}, RA::T, ZA::T, side::IMAS.BuildLayerSide) where {T<:Float64}

clip rails (rail_r, rail_z) so that they start/end along the intersection with the lines connecting
the magnetic axis (RA,ZA) and the point of maximum curvature along the equilibrium boundary (pr,pz)
"""
function clip_rails(rail_r::Vector{T}, rail_z::Vector{T}, pr::Vector{T}, pz::Vector{T}, RA::T, ZA::T, side::IMAS.BuildLayerSide) where {T<:Float64}
    IMAS.reorder_flux_surface!(rail_r, rail_z, RA, ZA)

    index = argmin(rail_r)
    rail_z = circshift(rail_z[1:end-1], index)
    rail_r = circshift(rail_r[1:end-1], index)

    # smooth out pr,pz to make sure curvature does not pick up anything spurious
    pr, pz = IMAS.resample_2d_path(pr, pz; method=:linear)
    pr, pz = IMAS.MXH(pr, pz, 2)(100)
    IMAS.reorder_flux_surface!(pr, pz, RA, ZA)

    theta = atan.(pz .- ZA, pr .- RA)
    weight = sin.(theta) .^ 2
    curve = abs.(IMAS.curvature(pr, pz)) .* weight

    index = Int[]
    index = pz .> ZA
    RU = @views pr[index][argmax(curve[index])]
    ZU = @views pz[index][argmax(curve[index])]
    index = pz .< ZA
    RL = @views pr[index][argmax(curve[index])]
    ZL = @views pz[index][argmax(curve[index])]

    α = 10.0

    ru = [RA, (RU - RA) * α + RA]
    zu = [ZA, (ZU - ZA) * α + ZA]
    intsc = IMAS.intersection(rail_r, rail_z, ru, zu)
    idx_u1 = intsc.indexes[1][1]
    crx_u1 = intsc.crossings[1]

    rl = [RA, (RL - RA) * α + RA]
    zl = [ZA, (ZL - ZA) * α + ZA]
    intsc = IMAS.intersection(rail_r, rail_z, rl, zl)
    idx_l2 = intsc.indexes[1][1]
    crx_l2 = intsc.crossings[1]

    if side in (_lfs_, _out_)
        rail_r = [crx_u1[1]; @views rail_r[idx_u1+1:idx_l2]; crx_l2[1]]
        rail_z = [crx_u1[2]; @views rail_z[idx_u1+1:idx_l2]; crx_l2[2]]

    elseif side == _hfs_
        for idx in idx_l2:-1:idx_u1
            popat!(rail_r, idx)
            popat!(rail_z, idx)
        end
        index = argmax(rail_r) - 1
        circshift!(rail_z, -index)
        circshift!(rail_r, -index)

    else
        error("clip_rails side can only be _hfs_ or _lfs_")
    end

    return rail_r, rail_z
end

"""
    init_pf_active!(
        pf_active::IMAS.pf_active,
        bd::IMAS.build,
        eqt::IMAS.equilibrium__time_slice;
        pf_coils_size::Union{Nothing,Float64,Vector{Float64}}=nothing,
        coils_cleareance::Union{Nothing,Float64,Vector{Float64}}=nothing)

Use build layers outline to initialize PF coils distribution
"""
function init_pf_active!(
    pf_active::IMAS.pf_active,
    bd::IMAS.build,
    eqt::IMAS.equilibrium__time_slice)

    OH_layer = IMAS.get_build_layer(bd.layer; type=_oh_)

    n_coils = [length(layer.coils_inside) for layer in bd.layer if !ismissing(layer, :coils_inside)]

    empty!(pf_active)
    resize!(bd.pf_active.rail, length(n_coils))

    krail = 1
    pf_coils_size = 0.0
    for (kl, layer) in enumerate(bd.layer)
        if ismissing(layer, :coils_inside) || layer.thickness == 0.0
            continue
        end

        if layer.type == Int(_oh_)

            # handle the case of OH inside of the TF
            # by creating a fake OH_layer that has a rectangular shape
            if OH_layer.side == Int(_hfs_)
                OH_layer = IMAS.freeze(OH_layer)
                height = 0.75 * maximum(OH_layer.outline.z) - minimum(OH_layer.outline.z)
                OH_layer.outline.r, OH_layer.outline.z = rectangle_shape(OH_layer.start_radius, OH_layer.end_radius, height)
                w_oh = OH_layer.thickness * 0.9
            else
                w_oh = maximum(OH_layer.outline.r) - minimum(OH_layer.outline.r)
            end
            coils_cleareance = w_oh / 8.0

            # OH coils are distributed on a rail within the OH region
            r_oh = sum(extrema(OH_layer.outline.r)) / 2.0

            z_ohcoils, h_oh = size_oh_coils(minimum(OH_layer.outline.z), maximum(OH_layer.outline.z), coils_cleareance[1], n_coils[1])
            bd.pf_active.rail[1].name = "OH"
            bd.pf_active.rail[1].coils_number = n_coils[1]
            bd.pf_active.rail[1].coils_cleareance = coils_cleareance[1]
            bd.pf_active.rail[1].outline.r = [r_oh, r_oh]
            bd.pf_active.rail[1].outline.z = [minimum(OH_layer.outline.z), maximum(OH_layer.outline.z)]
            bd.pf_active.rail[1].outline.distance = [-1.0, 1.0]
            for (kk, z_oh) in enumerate(z_ohcoils)
                k = length(pf_active.coil) + 1
                resize!(pf_active.coil, k)
                resize!(pf_active.coil[k].element, 1)
                coil = pf_active.coil[k]
                coil.identifier = "optim"
                coil.name = "OH $kk"
                coil.element[1].geometry.rectangle.r = r_oh
                coil.element[1].geometry.rectangle.z = z_oh
                coil.element[1].geometry.rectangle.width = w_oh
                coil.element[1].geometry.rectangle.height = h_oh
                coil.current.time = IMAS.top_ids(eqt).time
                coil.current.data = coil.current.time .* 0.0
                func = resize!(coil.function, :flux; wipe=false)
                func.description = "OH"
                func = resize!(coil.function, :shaping; wipe=false)
                func.description = "PF"
            end
            pf_coils_size = sqrt(w_oh * h_oh) * sqrt(2)

        else # PFs
            coils_cleareance = layer.thickness * 0.1
            pf_coils_size = layer.thickness * 0.9

            ngrid = 257
            if eqt.boundary.triangularity > 0.0
                dr = maximum(bd.layer[end].outline.r) / ngrid
            else
                rmask, zmask, mask = IMAS.structures_mask(bd; ngrid, layer_check=layer -> layer.material == "vacuum" || contains(lowercase(layer.name), "coils"))
                dr = (rmask[2] - rmask[1])
            end

            krail += 1
            nc = n_coils[krail]

            # add rail info to build IDS
            bd.pf_active.rail[krail].name = replace(replace(layer.name, "hfs " => ""), "lfs " => "")
            bd.pf_active.rail[krail].coils_number = nc
            bd.pf_active.rail[krail].coils_cleareance = coils_cleareance

            # generate rail between the two layers where coils will be placed and will be able to slide during the `optimization` phase
            if layer.side in (Int(_lfs_), Int(_out_))
                inner_layer = bd.layer[kl-1]
            else
                inner_layer = bd.layer[kl+1]
            end
            hull = IMAS.convex_hull(inner_layer.outline.r, inner_layer.outline.z; closed_polygon=true)
            rail0_r = [r for (r, z) in hull]
            rail0_z = [z for (r, z) in hull]
            if layer.side == Int(_out_)
                dcoil = (pf_coils_size + coils_cleareance) / 2
                rail0_r, rail0_z = buffer(rail0_r, rail0_z, dcoil)
            else
                layer_hfs = IMAS.get_build_layer(bd.layer; identifier=bd.layer[kl].identifier, fs=_hfs_)
                layer_lfs = IMAS.get_build_layer(bd.layer; identifier=bd.layer[kl].identifier, fs=_lfs_)
                dcoil = (layer_hfs.thickness + layer_lfs.thickness) / 2.0 / 2.0
                rail0_r, rail0_z = buffer(rail0_r, rail0_z, layer_hfs.thickness / 2.0, layer_lfs.thickness / 2.0)
            end
            coil_size = dcoil / sqrt(2) * 2
            rail0_r, rail0_z = IMAS.resample_2d_path(rail0_r, rail0_z; step=dr / 3.0)

            # let rails start along the lines connecting the magnetic axis and the point of maximum elongation
            RA = eqt.global_quantities.magnetic_axis.r
            ZA = eqt.global_quantities.magnetic_axis.z
            rail_r, rail_z = clip_rails(rail0_r, rail0_z, eqt.boundary.outline.r, eqt.boundary.outline.z, RA, ZA, IMAS.BuildLayerSide(layer.side))
            if layer.side == Int(_hfs_)
                index = rail_r .< maximum(rail_r) - 2.0 * dcoil
                valid_r, valid_z = rail_r[index], rail_z[index]
            else
                valid_r, valid_z = rail_r, rail_z
            end
            distance = cumsum(sqrt.(IMAS.gradient(valid_r) .^ 2 .+ IMAS.gradient(valid_z) .^ 2))

            # normalize distance between -1 and 1
            distance = (distance .- distance[1])
            distance = (distance ./ distance[end]) .* 2.0 .- 1.0

            # add rail info to build IDS
            bd.pf_active.rail[krail].outline.r = valid_r
            bd.pf_active.rail[krail].outline.z = valid_z
            bd.pf_active.rail[krail].outline.distance = distance

            if nc == 0
                continue
            end

            # uniformely distribute coils
            coils_distance = range(-1.0, 1.0, nc)
            r_coils = IMAS.interp1d(distance, valid_r).(coils_distance)
            z_coils = IMAS.interp1d(distance, valid_z).(coils_distance)
            z_coils = [abs(z) < 1E-6 ? 0.0 : z for z in z_coils]

            # populate IMAS data structure
            for (kk, (r, z)) in enumerate(zip(r_coils, z_coils))
                k = length(pf_active.coil) + 1
                resize!(pf_active.coil, k)
                resize!(pf_active.coil[k].element, 1)
                coil = pf_active.coil[k]
                coil.identifier = "optim"
                coil.name = "PF $kk"
                coil.element[1].geometry.rectangle.r = r
                coil.element[1].geometry.rectangle.z = z
                coil.element[1].geometry.rectangle.width = coil_size
                coil.element[1].geometry.rectangle.height = coil_size
                coil.current.time = IMAS.top_ids(eqt).time
                coil.current.data = coil.current.time .* 0.0
                bd.pf_active.rail[krail].name = "PF"
                func = resize!(coil.function, :shaping; wipe=false)
                func.description = "PF"
            end
        end
    end

    IMAS.set_coils_function(pf_active.coil, eqt.global_quantities.vacuum_toroidal_field.r0)

    return pf_active
end

function size_oh_coils(min_z::Real, max_z::Real, coils_cleareance::Real, coils_number::Int, height::Real=1.0, offset::Real=0.0)
    @assert 0.0 < height <= 1.0
    @assert -1.0 < offset < 1.0
    Δrail = max_z - min_z
    rail_offset = (max_z + min_z) / 2.0
    if coils_number == 1
        Δcoil = height * Δrail
        z = [rail_offset]
    else
        Δclear = coils_cleareance * coils_number
        Δcoil = (height * Δrail - Δclear) / coils_number
        z = range(-height * Δrail / 2.0 + Δcoil / 2.0, height * Δrail / 2.0 - Δcoil / 2.0, coils_number) .+ rail_offset
    end
    z = z .+ (offset * (1 - height) * Δrail)
    return z, Δcoil
end

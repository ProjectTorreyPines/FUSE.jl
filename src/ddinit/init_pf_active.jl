#= ================== =#
#  init pf_active IDS  #
#= ================== =#
"""
    init_pf_active!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd)

Initialize `dd.pf_active` starting from `ini` and `act` parameters
"""
function init_pf_active!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd)
    TimerOutputs.reset_timer!("init_pf_active")
    TimerOutputs.@timeit timer "init_pf_active" begin
        init_from = ini.general.init_from

        if init_from == :ods
            if length(dd1.pf_active.coil) > 0
                dd.pf_active = dd1.pf_active
                IMAS.set_coils_function(dd.pf_active.coil)

                # create rails from coils
                if false
                    OH = [coil for coil in dd.pf_active.coil if IMAS.is_ohmic_coil(coil)]
                    PF = [coil for coil in dd.pf_active.coil if !IMAS.is_ohmic_coil(coil)]
                    empty!(dd.build.pf_active.rail)
                    for (name, coils) in (("OH", OH), ("PF", PF))
                        rail = resize!(dd.build.pf_active.rail, length(dd.build.pf_active.rail) + 1)[end]
                        rail.name = name
                        rail.coils_cleareance = maximum([norm([coil.element[1].geometry.rectangle.width, coil.element[1].geometry.rectangle.height]) for coil in coils])
                        rail.coils_number = length(coils)
                        rail.outline.r = [coil.element[1].geometry.rectangle.r for coil in coils]
                        rail.outline.z = [coil.element[1].geometry.rectangle.z for coil in coils]
                        distance = cumsum(sqrt.(IMAS.gradient(rail.outline.r) .^ 2 .+ IMAS.gradient(rail.outline.z) .^ 2))
                        distance = (distance .- distance[1])
                        distance = (distance ./ distance[end]) .* 2.0 .- 1.0
                        rail.outline.distance = distance
                    end
                end

            else
                init_from = :scalars
            end
        end

        if init_from == :scalars
            n_coils = [ini.oh.n_coils]
            if any([contains(lowercase(layer.name), "coils") for layer in dd.build.layer])
                push!(n_coils, ini.pf_active.n_coils_inside)
            end
            if ini.pf_active.n_coils_outside > 0
                push!(n_coils, ini.pf_active.n_coils_outside)
            end
            init_pf_active!(dd.pf_active, dd.build, dd.equilibrium.time_slice[], n_coils)

        end

        coil_technology(dd.build.tf.technology, ini.tf.technology, :tf)
        coil_technology(dd.build.oh.technology, ini.oh.technology, :oh)
        coil_technology(dd.build.pf_active.technology, ini.pf_active.technology, :pf_active)

        return dd
    end
end

"""
    clip_rails(rail_r::Vector{T}, rail_z::Vector{T}, pr::Vector{T}, pz::Vector{T}, RA::T, ZA::T) where {T<:Real}

clip rails (rail_r, rail_z) so that they start/end along the intersection with the lines connecting
the magnetic axis (RA,ZA) and the point of maximum elongation along the equilibrium boundary (pr,pz)
"""
function clip_rails(rail_r::Vector{T}, rail_z::Vector{T}, pr::Vector{T}, pz::Vector{T}, RA::T, ZA::T) where {T<:Real}
    IMAS.reorder_flux_surface!(rail_r, rail_z)
    rail_z = circshift(rail_z[1:end-1], argmin(rail_r))
    rail_r = circshift(rail_r[1:end-1], argmin(rail_r))

    α = 1000.0
    index = []

    RU = pr[argmax(pz)]
    ZU = pz[argmax(pz)]
    ru = [RA, (RU - RA) * α + RA]
    zu = [ZA, (ZU - ZA) * α + ZA]
    push!(index, IMAS.intersection(rail_r, rail_z, ru, zu)[1][1][1])

    RL = pr[argmin(pz)]
    ZL = pz[argmin(pz)]
    rl = [RA, (RL - RA) * α + RA]
    zl = [ZA, (ZL - ZA) * α + ZA]
    push!(index, IMAS.intersection(rail_r, rail_z, rl, zl)[1][1][1])

    sort!(index)

    rail_r = rail_r[index[1]:index[2]]
    rail_z = rail_z[index[1]:index[2]]
    return rail_r, rail_z
end

"""
    init_pf_active!(
        pf_active::IMAS.pf_active,
        bd::IMAS.build,
        eqt::IMAS.equilibrium__time_slice,
        n_coils::Vector{TI};
        pf_coils_size::Union{Nothing,TR,Vector{TR}}=nothing,
        coils_cleareance::Union{Nothing,TR,Vector{TR}}=nothing) where {TI<:Int,TR<:Real}

Use build layers outline to initialize PF coils distribution
NOTE: n_coils

  - the first element in the array sets the number of coils in the OH
  - any subsequent element sets the number of coils for each of the vacuum regions in the build
"""
function init_pf_active!(
    pf_active::IMAS.pf_active,
    bd::IMAS.build,
    eqt::IMAS.equilibrium__time_slice,
    n_coils::Vector{TI};
    pf_coils_size::Union{Nothing,TR,Vector{TR}}=nothing,
    coils_cleareance::Union{Nothing,TR,Vector{TR}}=nothing) where {TI<:Int,TR<:Real}

    OH_layer = IMAS.get_build_layer(bd.layer; type=_oh_)

    empty!(pf_active)
    resize!(bd.pf_active.rail, length(n_coils))

    # coils_cleareance is an array the length of the rails
    if coils_cleareance === nothing
        coils_cleareance = (maximum(OH_layer.outline.r) - minimum(OH_layer.outline.r)) / 8.0
    end
    if isa(coils_cleareance, Number)
        coils_cleareance = [coils_cleareance for k in 1:length(n_coils)]
    end

    # OH coils are distributed on a rail within the OH region
    r_oh = sum(extrema(OH_layer.outline.r)) / 2.0
    w_oh = maximum(OH_layer.outline.r) - minimum(OH_layer.outline.r) - 2 * coils_cleareance[1]
    z_ohcoils, h_oh = size_oh_coils(OH_layer.outline.z, coils_cleareance[1], n_coils[1])
    bd.pf_active.rail[1].name = "OH"
    bd.pf_active.rail[1].coils_number = n_coils[1]
    bd.pf_active.rail[1].coils_cleareance = coils_cleareance[1]
    bd.pf_active.rail[1].outline.r = ones(length(z_ohcoils)) * r_oh
    bd.pf_active.rail[1].outline.z = z_ohcoils
    bd.pf_active.rail[1].outline.distance = range(-1, 1, n_coils[1])
    for z_oh in z_ohcoils
        k = length(pf_active.coil) + 1
        resize!(pf_active.coil, k)
        resize!(pf_active.coil[k].element, 1)
        pf_active.coil[k].identifier = "optim"
        pf_active.coil[k].name = "OH"
        pf_active.coil[k].element[1].geometry.rectangle.r = r_oh
        pf_active.coil[k].element[1].geometry.rectangle.z = z_oh
        pf_active.coil[k].element[1].geometry.rectangle.width = w_oh
        pf_active.coil[k].element[1].geometry.rectangle.height = h_oh
        @ddtime pf_active.coil[k].current.data = 0.0
    end

    # coils_cleareance is an array the length of the PF rails
    if pf_coils_size === nothing
        pf_coils_size = sqrt(w_oh * h_oh) * sqrt(2)
    end
    if isa(pf_coils_size, Number)
        pf_coils_size = [NaN, pf_coils_size, 1.5 * pf_coils_size]
    end

    # Now add actual PF coils to regions of vacuum
    gap_cryostat_index = [k for k in IMAS.get_build_indexes(bd.layer; fs=_out_) if bd.layer[k].material == "vacuum"][1]
    lfs_out_indexes = IMAS.get_build_indexes(bd.layer; fs=[_lfs_, _out_])
    krail = 1
    ngrid = 257
    rmask, zmask, mask = IMAS.structures_mask(bd; ngrid)
    dr = (rmask[2] - rmask[1])
    for k in lfs_out_indexes
        layer = bd.layer[k]

        if (k == gap_cryostat_index) && (length(n_coils) >= krail + 1) && (n_coils[krail+1] > 0)
            #pass
        elseif !contains(lowercase(layer.name), "coils")
            continue
        end

        krail += 1
        nc = n_coils[krail]

        # add rail info to build IDS
        bd.pf_active.rail[krail].name = replace(replace(layer.name, "hfs " => ""), "lfs " => "")
        bd.pf_active.rail[krail].coils_number = nc
        bd.pf_active.rail[krail].coils_cleareance = coils_cleareance[krail]

        # limit size of the pf_coils to fit in the vacuum region
        max_pf_inside_coil = (layer.end_radius - layer.start_radius - coils_cleareance[krail] * 2 * sqrt(2)) / sqrt(2)
        pf_coils_size[krail] = min(pf_coils_size[krail], max_pf_inside_coil)

        # generate rail between the two layers where coils will be placed and will be able to slide during the `optimization` phase
        coil_size = pf_coils_size[krail]
        dcoil = (coil_size + coils_cleareance[krail]) / 2 * sqrt(2)
        inner_layer = IMAS.get_build_layer(bd.layer; identifier=bd.layer[k-1].identifier, fs=_hfs_)
        rail_r, rail_z = buffer(inner_layer.outline.r, inner_layer.outline.z, dcoil)
        rail_r, rail_z = IMAS.resample_2d_path(rail_r, rail_z; step=dr / 3)

        # let rails start along the lines connecting the magnetic axis and the point of maximum elongation
        rail_r, rail_z = clip_rails(rail_r, rail_z, eqt.boundary.outline.r, eqt.boundary.outline.z, eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z)

        if true
            valid_r, valid_z = rail_r, rail_z
            distance = cumsum(sqrt.(IMAS.gradient(valid_r) .^ 2 .+ IMAS.gradient(valid_z) .^ 2))

        else
            # mark what regions on that rail do not intersect solid structures and can hold coils
            valid_k = []
            for (k, (r, z)) in enumerate(zip(rail_r, rail_z))
                ir = argmin(abs.(rmask .- r))
                iz = argmin(abs.(zmask .- z))
                if (ir < 1) || (ir > length(rmask)) || (iz < 1) || (iz > length(zmask))
                    continue
                end
                if (r > (minimum(rail_r) + (maximum(rail_r) - minimum(rail_r)) / 20)) && all(mask[ir, iz] .== 0)
                    push!(valid_k, k)
                end
            end
            if length(valid_k) == 0
                bd.pf_active.rail[krail].outline.r = Float64[]
                bd.pf_active.rail[krail].outline.z = Float64[]
                bd.pf_active.rail[krail].outline.distance = Float64[]
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
            valid_r = vcat(valid_r[istart+1:end], valid_r[1:istart])
            valid_z = vcat(valid_z[istart+1:end], valid_z[1:istart])

            # evaluate distance along rail
            d_distance = sqrt.(diff(vcat(valid_r, valid_r[1])) .^ 2.0 .+ diff(vcat(valid_z, valid_z[1])) .^ 2.0)
            d_distance[isnan.(d_distance)] .= 0.0
            valid_z = valid_z[d_distance.!=0]
            valid_r = valid_r[d_distance.!=0]
            distance = cumsum(d_distance)
            distance = distance[d_distance.!=0]
            # remove dcoil/sqrt(2) ends
            distance = (distance .- distance[1])
            valid_z = valid_z[distance.>dcoil/sqrt(2)]
            valid_r = valid_r[distance.>dcoil/sqrt(2)]
            distance = distance[distance.>dcoil/sqrt(2)]
            distance = (distance .- distance[end])
            valid_z = valid_z[distance.<-dcoil/sqrt(2)]
            valid_r = valid_r[distance.<-dcoil/sqrt(2)]
            distance = distance[distance.<-dcoil/sqrt(2)]
        end

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
            @ddtime pf_active.coil[k].current.data = 0.0
        end
    end

    IMAS.set_coils_function(pf_active.coil)

    return pf_active
end

function size_oh_coils(rail_outline_z, coils_cleareance, coils_number, height=1.0, offset=0.0)
    Δrail = maximum(rail_outline_z) - minimum(rail_outline_z)
    Δclear = coils_cleareance * coils_number
    Δcoil = (height * Δrail - Δclear) / coils_number
    rail_offset = (maximum(rail_outline_z) + minimum(rail_outline_z)) / 2.0
    z = range(-height * Δrail / 2.0 + Δcoil / 2.0, height * Δrail / 2.0 - Δcoil / 2.0, coils_number) .+ rail_offset
    z = z .+ (offset * (1 - height) * Δrail)
    return z, Δcoil
end

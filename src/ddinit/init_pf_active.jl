const coils_turns_spacing = 0.03

#= ================== =#
#  init pf_active IDS  #
#= ================== =#
function size_oh_coils(rail_outline_z, coils_cleareance, coils_number, height = 1.0, offset = 0.0)
    Δrail = maximum(rail_outline_z) - minimum(rail_outline_z)
    Δclear = coils_cleareance * coils_number
    Δcoil = (height * Δrail - Δclear) / coils_number
    rail_offset = (maximum(rail_outline_z) + minimum(rail_outline_z)) / 2.0
    z = LinRange(-height * Δrail / 2.0 + Δcoil / 2.0, height * Δrail / 2.0 - Δcoil / 2.0, coils_number) .+ rail_offset
    z = z .+ (offset * (1 - height) * Δrail)
    return z, Δcoil
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
function init_pf_active(
    pf_active::IMAS.pf_active,
    bd::IMAS.build,
    n_coils::Vector{TI};
    pf_coils_size::Union{Nothing,TR,Vector{TR}} = nothing,
    coils_cleareance::Union{Nothing,TR,Vector{TR}} = nothing,
    coils_elements_area::Union{Nothing,TR,Vector{TR}} = nothing) where {TI<:Int,TR<:Real}

    OH_layer = IMAS.get_build(bd, type = 1)

    empty!(pf_active)
    resize!(bd.pf_active.rail, length(n_coils))

    # coils_cleareance is an array the lenght of the rails
    if coils_elements_area === nothing
        coils_elements_area = 0.0025
    end
    if isa(coils_elements_area, Number)
        coils_elements_area = [coils_elements_area for k in 1:length(n_coils)]
    end

    # coils_cleareance is an array the lenght of the rails
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
    bd.pf_active.rail[1].coils_elements_area = coils_elements_area[1]
    bd.pf_active.rail[1].coils_cleareance = coils_cleareance[1]
    bd.pf_active.rail[1].outline.r = ones(length(z_ohcoils)) * r_oh
    bd.pf_active.rail[1].outline.z = z_ohcoils
    bd.pf_active.rail[1].outline.distance = range(-1, 1, length = n_coils[1])
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
        set_turns_from_spacing!(pf_active.coil[k], coils_turns_spacing, +1)
        @ddtime pf_active.coil[k].current.data = 0.0
    end

    # coils_cleareance is an array the lenght of the PF rails
    if pf_coils_size === nothing
        pf_coils_size = sqrt(w_oh * h_oh)
    end
    if isa(pf_coils_size, Number)
        pf_coils_size = [NaN, pf_coils_size, 1.5 * pf_coils_size]
    end

    # Now add actual PF coils to regions of vacuum
    krail = 1
    ngrid = 257
    rmask, zmask, mask = IMAS.structures_mask(bd, ngrid = ngrid)
    dr = (rmask[2] - rmask[1])
    for (k, layer) in enumerate(bd.layer)
        if ((layer.hfs == -1) || (k == length(bd.layer))) && !ismissing(layer, :material) && (lowercase(layer.material) == "vacuum")

            krail += 1
            nc = n_coils[krail]

            # add rail info to build IDS
            bd.pf_active.rail[krail].name = replace(replace(layer.name, "hfs " => ""), "lfs " => "")
            bd.pf_active.rail[krail].coils_number = nc
            bd.pf_active.rail[krail].coils_elements_area = coils_elements_area[krail]
            bd.pf_active.rail[krail].coils_cleareance = coils_cleareance[krail]

            # limit size of the pf_coils to fit in the vacuum region
            max_pf_inside_coil = (layer.end_radius - layer.start_radius - coils_cleareance[krail] * 2 * sqrt(2)) / sqrt(2)
            pf_coils_size[krail] = min(pf_coils_size[krail], max_pf_inside_coil)

            # generate rail between the two layers where coils will be placed and will be able to slide during the `optimization` phase
            coil_size = pf_coils_size[krail]
            dcoil = (coil_size + coils_cleareance[krail]) / 2 * sqrt(2)
            inner_layer = IMAS.get_build(bd, identifier = bd.layer[k-1].identifier, hfs = 1)
            poly = LibGEOS.buffer(xy_polygon(inner_layer.outline.r, inner_layer.outline.z), dcoil)
            rail_r = [v[1] for v in LibGEOS.coordinates(poly)[1]]
            rail_z = [v[2] for v in LibGEOS.coordinates(poly)[1]]
            rail_r, rail_z = IMAS.resample_2d_line(rail_r, rail_z, dr / 3)

            # mark what regions on that rail do not intersect solid structures and can hold coils
            valid_k = []
            for (k, (r, z)) in enumerate(zip(rail_r, rail_z))
                ir = argmin(abs.(rmask .- r))
                iz = argmin(abs.(zmask .- z))
                if (ir < 1) || (ir > length(rmask)) || (iz < 1) || (iz > length(zmask))
                    continue
                end
                if all(mask[ir, iz] .== 0)
                    push!(valid_k, k)
                end
            end
            if length(valid_k) == 0
                bd.pf_active.rail[krail].outline.r = Real[]
                bd.pf_active.rail[krail].outline.z = Real[]
                bd.pf_active.rail[krail].outline.distance = Real[]
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
            coils_distance = range(-(1 - 0.25 / nc), 1 - 0.25 / nc, length = nc)
            r_coils = IMAS.interp1d(distance, valid_r).(coils_distance)
            z_coils = IMAS.interp1d(distance, valid_z).(coils_distance)
            z_coils = [abs(z) < 1E-6 ? 0 : z for z in z_coils]

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

    return pf_active
end

function init_pf_active(dd::IMAS.dd, par::Parameters)
    init_from = par.general.init_from

    if init_from == :gasc
        init_from = :scalars

    elseif init_from == :ods
        dd1 = IMAS.json2imas(par.ods.filename)
        if length(keys(dd1.pf_active)) > 0
            dd.global_time = max(dd.global_time, maximum(dd1.pf_active.time))
            dd.pf_active = dd1.pf_active
        else
            init_from = :scalars
        end
    end

    if init_from == :scalars
        n_coils = [par.pf_active.n_oh_coils]
        if par.pf_active.n_pf_coils_inside > 0
            push!(n_coils, par.pf_active.n_pf_coils_inside)
        end
        push!(n_coils, par.pf_active.n_pf_coils_outside)
        init_pf_active(dd.pf_active, dd.build, n_coils)
    end

    assign_coils_materials(dd, par)

    return dd
end

function assign_coils_materials(dd::IMAS.dd, par::Parameters)
    for coil_type in [:tf, :oh, :pf_active]
        coil_tech = getproperty(dd.build, coil_type).technology
        for property in fieldnames(IMAS.build__tf__technology)
            if property == :_parent
                continue
            end
            coil_params = getproperty(par, coil_type).technology
            if !ismissing(coil_params, property)
                setproperty!(coil_tech, property, getproperty(coil_params, property))
            end
        end
    end
end

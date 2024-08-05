#= ========== =#
#  PF passive  #
#= ========== =#
Base.@kwdef mutable struct FUSEparameters__ActorPassiveStructures{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    do_plot::Entry{Bool} = act_common_parameters(do_plot=false)
end

mutable struct ActorPassiveStructures{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorPassiveStructures{P}
    function ActorPassiveStructures(dd::IMAS.dd{D}, par::FUSEparameters__ActorPassiveStructures{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorPassiveStructures)
        par = par(kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorPassiveStructures(dd::IMAS.dd, act::ParametersAllActors; kw...)

Populates `pf_passive` structures
"""
function ActorPassiveStructures(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorPassiveStructures(dd, act.ActorPassiveStructures; kw...)
    step(actor)
    finalize(actor)
    if actor.par.do_plot
        display(plot(dd.pf_passive))
    end
    return actor
end

function _step(actor::ActorPassiveStructures)
    dd = actor.dd

    empty!(dd.pf_passive)

    # all LFS layers that do not intersect with structures
    ilayers = IMAS.get_build_indexes(dd.build.layer; fs=IMAS._lfs_)
    ilayers = vcat(ilayers[1] - 1, ilayers)
    for k in ilayers
        l = dd.build.layer[k]
        l1 = dd.build.layer[k+1]
        if l1.material == "vacuum"
            continue
        end
        if all(l1.type != structure.type for structure in dd.build.structure)
            poly = IMAS.join_outlines(l.outline.r, l.outline.z, l1.outline.r, l1.outline.z)
            add_pf_passive_loop(dd.pf_passive, l1, poly[1], poly[2])
        end
    end

    # then all structures
    for structure in dd.build.structure
        add_pf_passive_loop(dd.pf_passive, structure)
    end

    # OH and plug
    add_pf_passive_loop(dd.pf_passive, dd.build.layer[2])
    if dd.build.layer[1].material != "vacuum"
        add_pf_passive_loop(dd.pf_passive, dd.build.layer[1])
    end

    # cryostat
    layer = dd.build.layer[end]
    if layer.type == Int(IMAS._cryostat_)
        i = argmin(abs.(layer.outline.r))
        r = vcat(layer.outline.r[i+1:end-1], layer.outline.r[1:i])
        z = vcat(layer.outline.z[i+1:end-1], layer.outline.z[1:i])
        layer = dd.build.layer[end-1]
        i = argmin(abs.(layer.outline.r))
        r = vcat(layer.outline.r[i+1:end-1], layer.outline.r[1:i], reverse(r), layer.outline.r[i+1:i+1])
        z = vcat(layer.outline.z[i+1:end-1], layer.outline.z[1:i], reverse(z), layer.outline.z[i+1:i+1])
        add_pf_passive_loop(dd.pf_passive, dd.build.layer[end], r, z)
    end

    return actor
end

function add_pf_passive_loop(pf_passive::IMAS.pf_passive, layer::IMAS.build__layer)
    return add_pf_passive_loop(pf_passive, layer, layer.outline.r, layer.outline.z)
end

function add_pf_passive_loop(pf_passive::IMAS.pf_passive, layer::IMAS.build__layer, r::Vector{T}, z::Vector{T}) where {T<:Real}
    return add_pf_passive_loop(pf_passive, layer.name, "layer $(IMAS.index(layer))", r, z)
end

function add_pf_passive_loop(pf_passive::IMAS.pf_passive, structure::IMAS.build__structure)
    return add_pf_passive_loop(pf_passive, structure.name, "structure $(IMAS.index(structure))", structure.outline.r, structure.outline.z)
end

function add_pf_passive_loop(pf_passive::IMAS.pf_passive, name::AbstractString, identifier::AbstractString, r::AbstractVector{T}, z::AbstractVector{T}) where {T<:Real}
    resize!(resize!(pf_passive.loop, length(pf_passive.loop) + 1)[end].element, 1)
    pf_passive.loop[end].name = replace(name, r"[lh]fs " => "")
    pf_passive.loop[end].element[end].geometry.outline.r = r
    pf_passive.loop[end].element[end].geometry.outline.z = z
    pf_passive.loop[end].element[end].geometry.geometry_type = IMAS.name_2_index(pf_passive.loop[end].element[end].geometry)[:outline]
    if !isempty(identifier)
        pf_passive.loop[end].element[end].identifier = identifier
    end
    return pf_passive.loop[end]
end

"""
    layer_quads(inner_layer::IMAS.build__layer, outer_layer::IMAS.build__layer, precision::Float64, min_n_segments::Int)

Build quads between two layers
"""
function layer_quads(inner_layer::IMAS.build__layer, outer_layer::IMAS.build__layer, precision::Float64, min_n_segments::Int)
    inner_outline = IMAS.closed_polygon(inner_layer.outline.r, inner_layer.outline.z)
    outer_outline = IMAS.closed_polygon(outer_layer.outline.r, outer_layer.outline.z)

    # reorder surface so that it starts on the hfs
    pr = inner_outline.r
    pz = inner_outline.z
    R0 = (maximum(pr) + minimum(pr)) * 0.5
    Z0 = (maximum(pz) + minimum(pz)) * 0.5
    indexes, crossings = IMAS.intersection(pr, pz, [0.0, R0], [Z0, Z0])
    @views pr = [pr[1:indexes[1][1]]; crossings[1][1]; pr[indexes[1][1]+1:end]]
    @views pz = [pz[1:indexes[1][1]]; crossings[1][2]; pz[indexes[1][1]+1:end]]
    istart = indexes[1][1] + 1
    IMAS.reorder_flux_surface!(pr, pz, istart)

    # simplify surface polygon
    @views R1, Z1 = IMAS.rdp_simplify_2d_path(pr[1:end-1], pz[1:end-1], precision)
    @views R1 = R1[2:end]
    @views Z1 = Z1[2:end]

    # split long segments
    L_inner = sum(sqrt.(diff(inner_layer.outline.r).^2 .+ diff(inner_layer.outline.z).^2))
    max_seg_length = L_inner / min_n_segments
    R1, Z1 = IMAS.split_long_segments(R1, Z1, max_seg_length)

    # radiate lines from polygon vertices
    rays = IMAS.polygon_rays(collect(zip(R1, Z1)), -100.0, 100.0)

    # generate vertices of inner and outer surfaces
    qR1 = Float64[]
    qZ1 = Float64[]
    qR2 = Float64[]
    qZ2 = Float64[]
    for (k, ray) in enumerate(rays)
        _, inner_crossings = IMAS.intersection(inner_outline.r, inner_outline.z, [ray[1][1], ray[2][1]], [ray[1][2], ray[2][2]])
        _, outer_crossings = IMAS.intersection(outer_outline.r, outer_outline.z, [ray[1][1], ray[2][1]], [ray[1][2], ray[2][2]])
        if isempty(inner_crossings) || isempty(outer_crossings)
            continue
        end
        d = [sqrt.((c[1] - R1[k])^2 + (c[2] - Z1[k])^2) for c in inner_crossings]
        i = argmin(d)
        push!(qR1, inner_crossings[i][1])
        push!(qZ1, inner_crossings[i][2])
        d = [sqrt.((c[1] - R1[k])^2 + (c[2] - Z1[k])^2) for c in outer_crossings]
        i = argmin(d)
        push!(qR2, outer_crossings[i][1])
        push!(qZ2, outer_crossings[i][2])
    end

    # define quads
    quads = [define_quad(ka, qR1, qR2, qZ1, qZ2) for ka in eachindex(qR1)]

    return quads
end

function define_quad(ka, qR1, qR2, qZ1, qZ2)
    kb = IMAS.getindex_circular(1:length(qR1), ka + 1)
    rr = @SVector[qR1[ka], qR1[kb], qR2[kb], qR2[ka]]
    zz = @SVector[qZ1[ka], qZ1[kb], qZ2[kb], qZ2[ka]]
    RC = sum(rr) * 0.25
    ZC = sum(zz) * 0.25
    index = sortperm(atan.(zz .- ZC, rr .- RC))
    return (rr[index], zz[index])
end
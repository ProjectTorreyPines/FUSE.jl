#= ========== =#
#  PF passive  #
#= ========== =#
Base.@kwdef mutable struct FUSEparameters__ActorPassiveStructures{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    #== actor parameters ==#
    wall_precision::Entry{Float64} = Entry{Float64}("-", "Precision for making wall quadralaterals"; default=1.0)
    min_n_segments::Entry{Int} = Entry{Int}("-", "Minimum number of quadralaterals"; default=15)
    #== display and debugging parameters ==#
    do_plot::Entry{Bool} = act_common_parameters(do_plot=false)
end

mutable struct ActorPassiveStructures{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorPassiveStructures{P}}
    function ActorPassiveStructures(dd::IMAS.dd{D}, par::FUSEparameters__ActorPassiveStructures{P}; kw...) where {D<:Real,P<:Real}
        logging_actor_init(ActorPassiveStructures)
        par = OverrideParameters(par; kw...)
        return new{D,P}(dd, par)
    end
end

"""
    ActorPassiveStructures(dd::IMAS.dd, act::ParametersAllActors; kw...)

Populates `pf_passive` IDS based on the vacuum vessel layer(s)
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
    par = actor.par

    empty!(dd.pf_passive)

    # The vacuum vessel can have multiple layers
    kvessels = IMAS.get_build_indexes(dd.build.layer; type=_vessel_, fs=_lfs_)
    if isempty(kvessels)
        @warn "No vessel found. Can't compute vertical stability metrics"
        return actor
    end
    for kvessel in kvessels
        kout = kvessel
        kin = kvessel - 1

        # resistivity just takes the material from the outermost build layer;
        # It does not account for toroidal breaks, heterogeneous materials,
        # or builds with "water" vacuum vessels
        mat_vv = Material(dd.build.layer[kout].material)
        if ismissing(mat_vv) || ismissing(mat_vv.electrical_conductivity)
            mat_vv = Material(:steel)
        end
        resistivity = 1.0 / mat_vv.electrical_conductivity(; temperature=273.15)

        thickness_to_radius = dd.build.layer[kin].thickness / (2π * (dd.build.layer[kin].start_radius - IMAS.opposite_side_layer(dd.build.layer[kin]).end_radius) / 2)

        quads = layer_quads(dd.build.layer[kin], dd.build.layer[kout], par.wall_precision * thickness_to_radius, par.min_n_segments)
        for (k,quad) in enumerate(quads)
            add_pf_passive_loop(dd.pf_passive, dd.build.layer[kout].name, "VV $k", quad[1], quad[2]; resistivity)
        end
    end

    return actor
end

function add_pf_passive_loop(pf_passive::IMAS.pf_passive, name::AbstractString, identifier::AbstractString, r::AbstractVector{T}, z::AbstractVector{T}; resistivity::T) where {T<:Real}
    resize!(resize!(pf_passive.loop, length(pf_passive.loop) + 1)[end].element, 1)
    loop = pf_passive.loop[end]
    loop.name = replace(name, r"[lh]fs " => "")
    loop.element[end].geometry.outline.r = r
    loop.element[end].geometry.outline.z = z
    loop.element[end].geometry.geometry_type = IMAS.name_2_index(loop.element[end].geometry)[:outline]
    if !isempty(identifier)
        loop.element[end].identifier = identifier
    end
    loop.resistivity = resistivity
    return loop
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
    L_inner = IMAS.perimeter(inner_layer.outline.r, inner_layer.outline.z)
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
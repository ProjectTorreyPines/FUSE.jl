#= ========== =#
#  PF passive  #
#= ========== =#
Base.@kwdef mutable struct FUSEparameters__ActorPassiveStructures{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    do_plot::Entry{Bool} = Entry{Bool}("-", "Plot"; default=false)
end

mutable struct ActorPassiveStructures{D,P} <: ReactorAbstractActor
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
    ilayers = IMAS.get_build_indexes(dd.build.layer, fs=IMAS._lfs_)
    ilayers = vcat(ilayers[1] - 1, ilayers)
    for k in ilayers
        l = dd.build.layer[k]
        l1 = dd.build.layer[k+1]
        if l1.material == "Vacuum"
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
    if dd.build.layer[1].material != "Vacuum"
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
    if !isempty(identifier)
        pf_passive.loop[end].element[end].identifier = identifier
    end
    return pf_passive.loop[end]
end
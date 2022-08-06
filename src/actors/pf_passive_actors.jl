#= ========== =#
#  PF passive  #
#= ========== =#

mutable struct ActorPassiveStructures <: ReactorAbstractActor
    dd::IMAS.dd
    par::ParametersActor
end

function ParametersActor(::Type{Val{:ActorPassiveStructures}})
    par = ParametersActor(nothing)
    par.do_plot = Entry(Bool, "", "plot"; default=false)
    return par
end

"""
    ActorPassiveStructures(dd::IMAS.dd, act::ParametersAllActors; kw...)

Populates `pf_passive` structures
"""
function ActorPassiveStructures(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorPassiveStructures(kw...)
    actor = ActorPassiveStructures(dd, par)
    step(actor)
    finalize(actor)
    if par.do_plot
        display(plot(dd.pf_passive))
    end
    return actor
end

function ActorPassiveStructures(dd::IMAS.dd, par::ParametersActor; kw...)
    par = par(kw...)
    ActorPassiveStructures(dd, par)
end

function step(actor::ActorPassiveStructures)
    dd = actor.dd

    # all LFS
    ilayers = IMAS.get_build(dd.build, fs=IMAS._lfs_, return_only_one=false, return_index=true)
    ilayers = vcat(ilayers[1] - 1, ilayers)
    for k in ilayers
        l = dd.build.layer[k]
        l1 = dd.build.layer[k+1]
        if l1.material == "Vacuum"
            continue
        end
        poly = IMAS.join_outlines(l.outline.r, l.outline.z, l1.outline.r, l1.outline.z)
        add_pf_passive_loop(dd.pf_passive, l1, poly[1], poly[2])
    end

    # OH and plug
    add_pf_passive_loop(dd.pf_passive, dd.build.layer[2], dd.build.layer[2].outline.r, dd.build.layer[2].outline.z)
    if dd.build.layer[1].material != "Vacuum"
        add_pf_passive_loop(dd.pf_passive, dd.build.layer[1], dd.build.layer[1].outline.r, dd.build.layer[1].outline.z)
    end

    # cryostat
    layer = dd.build.layer[end]
    if layer.type == Int(IMAS._cryostat_)
        i = findfirst(layer.outline.r .== 0)
        r = vcat(layer.outline.r[i+1:end-1], layer.outline.r[1:i])
        z = vcat(layer.outline.z[i+1:end-1], layer.outline.z[1:i])
        layer = dd.build.layer[end-1]
        i = findfirst(layer.outline.r .== 0)
        r = vcat(layer.outline.r[i+1:end-1], layer.outline.r[1:i], reverse(r), layer.outline.r[i+1:i+1])
        z = vcat(layer.outline.z[i+1:end-1], layer.outline.z[1:i], reverse(z), layer.outline.z[i+1:i+1])
        add_pf_passive_loop(dd.pf_passive, dd.build.layer[end], r, z)
    end
end

function add_pf_passive_loop(pf_passive::IMAS.pf_passive, layer::IMAS.build__layer, r::AbstractVector{T}, z::AbstractVector{T}) where {T<:Real}
    resize!(resize!(pf_passive.loop, length(pf_passive.loop) + 1)[end].element, 1)
    pf_passive.loop[end].name = replace(layer.name, r"[lh]fs " => "")
    pf_passive.loop[end].element[end].geometry.outline.r = r
    pf_passive.loop[end].element[end].geometry.outline.z = z
end
#= ============== =#
#  materials cost #
#= ============== =#
#NOTE: material should be priced by Kg
#NOTE: if something is priced by m^3 then it is for a specific part already
function unit_cost(material::String)
    if material == "Vacuum"
        return 0.0 # $M/m^3
    elseif material == "ReBCO"
        return 87.5/2 # $M/m^3
    elseif material == "Nb3Sn"
        return 1.66 # $M/m^3
    elseif contains(lowercase(material), "steel")
        return 0.36 # $M/m^3
    elseif material == "Tungsten"
        return 0.36 # $M/m^3
    elseif material == "Copper"
        return 0.5 # $M/m^3
    elseif material == "Water, Liquid"
        return 0.0 # $M/m^3
    elseif material == "lithium-lead"
        return 0.75 # $M/m^3
    elseif material == "FLiBe"
        return 0.75*3 # $M/m^3
    elseif contains(lowercase(material), "plasma")
        return 0.0 # $M/m^3
    else
        error("Material `$material` has no price \$M/m³")
    end
end

function unit_cost(coil_tech::Union{IMAS.build__tf__technology,IMAS.build__oh__technology,IMAS.build__pf_active__technology})
    if coil_tech.material == "Copper"
        return unit_cost("Copper")
    else
        fraction_cable = 1 - coil_tech.fraction_stainless - coil_tech.fraction_void
        fraction_SC = fraction_cable * coil_tech.ratio_SC_to_copper
        fraction_copper = fraction_cable - fraction_SC
        (coil_tech.fraction_stainless * unit_cost("Steel, Stainless 316") + fraction_copper * unit_cost("Copper") + fraction_SC * unit_cost(coil_tech.material))
    end
end


#= ======= =#
#  costing #
#= ======= =#

mutable struct CostingActor <: AbstractActor
    dd::IMAS.dd
end

function ActorParameters(::Type{Val{:CostingActor}})
    par = ActorParameters(nothing)
    return par
end

function CostingActor(dd::IMAS.dd, act::ActorParameters; kw...)
    par = act.CostingActor(kw...)
    actor = CostingActor(dd)
    step(actor)
    finalize(actor)
    return actor
end

function step(actor::CostingActor)
    return actor
end

function cost(layer::IMAS.build__layer)
    if layer.type == Int(_oh_)
        build = IMAS.parent(IMAS.parent(layer))
        unit_cost(build.oh.technology) * layer.volume
    elseif layer.type == Int(_tf_)
        build = IMAS.parent(IMAS.parent(layer))
        unit_cost(build.tf.technology) * layer.volume
    else
        unit_cost(layer.material) * layer.volume
    end
end
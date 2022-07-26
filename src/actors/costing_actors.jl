#= ============== =#
#  materials cost #
#= ============== =#
#NOTE: material should be priced by Kg
#NOTE: if something is priced by m^3 then it is for a specific part already
function unit_cost(material::String)
    if material == "Vacuum"
        return 0.0 # $M/m^3
    elseif material == "ReBCO"
        return 87.5 / 2 # $M/m^3
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
        return 0.75 * 3 # $M/m^3
    elseif contains(lowercase(material), "plasma")
        return 0.0 # $M/m^3
    else
        error("Material `$material` has no price \$M/m³")
    end
end

function cost(layer::IMAS.build__layer)
    if layer.type == Int(_oh_)
        build = IMAS.parent(IMAS.parent(layer))
        return unit_cost(build.oh.technology) * layer.volume
    elseif layer.type == Int(_tf_)
        build = IMAS.parent(IMAS.parent(layer))
        return unit_cost(build.tf.technology) * layer.volume
    elseif layer.type == Int(_shield_)
        return layer.volume * 0.29  # $M/m^3
    elseif layer.type == Int(_blanket_)
        return layer.volume * 0.75  # $M/m^3
    elseif layer.type ∈ [Int(_wall_), Int(_vessel_), Int(_cryostat_)]
        return layer.volume * 0.36  # $M/m^3
    else
        return unit_cost(layer.material) * layer.volume
    end
end

function cost(ecl::IMAS.ec_launchers__launcher)
    ecl.available_launch_power / 1E6 * 3.0 # $/W #ARIES
end

function cost(ica::IMAS.ic_antennas__antenna)
    ica.available_launch_power / 1E6 * 1.64 #$/W ARIES
end

function cost(lha::IMAS.lh_antennas__antenna)
    lha.available_launch_power / 1E6 * 2.13 #$/W ARIES
end

function cost(nbu::IMAS.nbi__unit)
    nbu.available_launch_power / 1E6 * 4.93 #$/W ARIES
end

function unit_cost(coil_tech::Union{IMAS.build__tf__technology,IMAS.build__oh__technology,IMAS.build__pf_active__technology})
    if coil_tech.material == "Copper"
        return unit_cost("Copper")
    else
        fraction_cable = 1 - coil_tech.fraction_stainless - coil_tech.fraction_void
        fraction_SC = fraction_cable * coil_tech.ratio_SC_to_copper
        fraction_copper = fraction_cable - fraction_SC
        return (coil_tech.fraction_stainless * unit_cost("Steel, Stainless 316") + fraction_copper * unit_cost("Copper") + fraction_SC * unit_cost(coil_tech.material))
    end
end

function cost(pf_active::IMAS.pf_active)
    dd = IMAS.top_dd(pf_active)
    c = Dict("OH" => 0.0, "PF" => 0.0)
    for coil in pf_active.coil
        if coil.name == "OH"
            c["OH"] += cost(coil, dd.build.oh.technology)
        else
            c["PF"] += cost(coil, dd.build.pf_active.technology)
        end
    end
    return c
end

function cost(coil::IMAS.pf_active__coil, technology::Union{IMAS.build__tf__technology,IMAS.build__oh__technology,IMAS.build__pf_active__technology})
    return IMAS.volume(coil) * unit_cost(technology)
end

#= ============ =#
#  ActorCosting  #
#= ============ =#

Base.@kwdef mutable struct ActorCosting <: FacilityAbstractActor
    dd::IMAS.dd
    par::ParametersActor
end

function ParametersActor(::Type{Val{:ActorCosting}})
    par = ParametersActor(nothing)
    return par
end

"""
    ActorCosting(dd::IMAS.dd, act::ParametersAllActors; kw...)

This actor estimates the cost of the fusion power plant.

!!! note 
    Stores data in `dd.costing`
"""
function ActorCosting(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorCosting(kw...)
    actor = ActorCosting(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorCosting(dd::IMAS.dd, par::ParametersActor; kw...)
    par = par(kw...)
    return ActorCosting(dd, par)
end

function step(actor::ActorCosting)
    dd = actor.dd
    cst = dd.costing
    empty!(cst)

    # TOKAMAK
    sys = resize!(cst.system, "name" => "tokamak")
    for layer in dd.build.layer
        if layer.fs == Int(_lfs_)
            continue # avoid double counting of hfs and lfs layers
        elseif layer.type == Int(_oh_)
            continue # avoid double counting of oh
        end
        c = cost(layer)
        if c > 0
            sub = resize!(sys.subsystem, "name" => replace(layer.name, r"^hfs " => ""))
            sub.cost = c
        end
    end
    for (name, c) in cost(dd.pf_active)
        sub = resize!(sys.subsystem, "name" => name)
        sub.cost = c
    end

    # HCD
    sys = resize!(cst.system, "name" => "hcd")
    for hcd in vcat(dd.ec_launchers.launcher, dd.ic_antennas.antenna, dd.lh_antennas.antenna, dd.nbi.unit)
        c = cost(hcd)
        if c > 0
            sub = resize!(sys.subsystem, "name" => hcd.name)
            sub.cost = c
        end
    end

    return actor
end

function finalize(actor::ActorCosting)
    # sort system/subsystem costs
    sort!(actor.dd.costing.system, by=x -> x.cost, rev=true)
    for sys in actor.dd.costing.system
        sort!(sys.subsystem, by=x -> x.cost, rev=true)
    end
end
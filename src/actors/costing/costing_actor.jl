function today_dollars(dollars::Real, year::Integer, generation::Symbol, costing::IMAS.costing)
    # table of inflation from dawn of time till doday + dd.costing.future_inflation_rate for extrapolation
    return dollars
end

include("ARIES_costing.jl")

include("Sheffield_costing.jl")

include("costing.jl")

#= ============ =#
#  ActorCosting  #
#= ============ =#
Base.@kwdef mutable struct FUSEparameters__ActorCosting{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    model::Switch{Symbol} = Switch(Symbol, [:FUSE, :ARIES, :Sheffield], "-", "Costing model"; default=:FUSE)
    construction_start_date::
    construction_lead_time::
    dd.costing.generation::Switch(Symbol, [:prototype, :1stofkind, :nthofkind], ...)
    future_inflation_rate::
    land_space::Entry{T} = Entry(T, "acres", "Plant site space required in acres"; default=1000.0)
    building_volume::Entry{T} = Entry(T, "m^3", "Volume of the tokmak building"; default=140.0e3)
    interest_rate::Entry{T} = Entry(T, "-", "Anual interest rate fraction of direct capital cost"; default=0.05)
    indirect_cost_rate::Entry{T} = Entry(T, "-", "Indirect cost associated with construction, equipment, services, energineering construction management and owners cost"; default=0.4)
    lifetime::Entry{Int} = Entry(Int, "years", "lifetime of the plant"; default=40)
    availability::Entry{T} = Entry(T, "-", "availability fraction of the plant"; default=0.803)
    escalation_fraction::Entry{T} = Entry(T, "-", "yearly escalation fraction based on risk assessment"; default=0.05)
    blanket_lifetime::Entry{T} = Entry(T, "years", "lifetime of the blanket"; default=6.8)
end

mutable struct ActorCosting <: FacilityAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorCosting
    function ActorCosting(dd::IMAS.dd, par::FUSEparameters__ActorCosting; kw...)
        logging_actor_init(ActorCosting)
        par = par(kw...)
        return new(dd, par)
    end
end

"""
    ActorCosting(dd::IMAS.dd, act::ParametersAllActors; kw...)

Estimates the cost of building, operating, and recommission the fusion power plant.

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

function _step(actor::ActorCosting)
    par = actor.par
    dd = actor.dd
    cst = dd.costing
    cost_direct = cst.cost_direct_capital
    cost_ops = cst.cost_operations
    cost_decom = cst.cost_decommissioning

    empty!(cst)

    dd.costing.timeline.construction.start_date = par.construction_start_date
    dd.costing.timeline.construction.lead_time = par.construction_lead_time
    dd.costing.timeline.lifetime = par.lifetime
    dd.costing.model = 
    dd.costing.generation = 
    dd.costing.future_inflation_rate =

    if par.model == :ARIES
        costing_ARIES(dd, par)
    elseif par.model == :Sheffield
        costing_Sheffield(dd, par)
    elseif par.model == :FUSE
        costing_FUSE(dd, par)
    end

    ###### Levelized Cost Of Electricity  ######
    capital_cost_rate = par.interest_rate / (1 - (1 + par.interest_rate)^(-1.0 * par.lifetime))
    lifetime_cost = 0.0
    for year in 1:par.lifetime
        yearly_cost = (capital_cost_rate * cost_direct.cost + cost_ops.yearly_cost + cost_decom.cost / par.lifetime)
        lifetime_cost += (1.0 + par.escalation_fraction) * (1.0 + par.indirect_cost_rate) * yearly_cost
    end
    dd.costing.cost_lifetime = lifetime_cost
    dd.costing.levelized_CoE = (dd.costing.cost_lifetime * 1E6) / (par.lifetime * 24 * 365 * power_electric_net / 1e3 * par.availability)
    return actor
end

function _finalize(actor::ActorCosting)
    # sort system/subsystem by their costs
    sort!(actor.dd.costing.cost_direct_capital.system, by=x -> x.cost, rev=true)
    for sys in actor.dd.costing.cost_direct_capital.system
        sort!(sys.subsystem, by=x -> x.cost, rev=true)
    end
end

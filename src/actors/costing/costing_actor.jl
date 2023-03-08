function today_dollars(dollars::Real, year::Integer, generation::Symbol, costing::IMAS.costing)
    # table of inflation from dawn of time till doday + dd.costing.future_inflation_rate for extrapolation
    return dollars
end 

#= ============ =#
#  ActorCosting  #
#= ============ =#
Base.@kwdef mutable struct FUSEparameters__ActorCosting{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    model::Switch{Symbol} = Switch(Symbol, [:FUSE, :ARIES, :Sheffield], "-", "Costing model"; default=:FUSE)
    #construction_start_date::
    #construction_lead_time::
    #dd.costing.generation::Switch(Symbol, [:prototype, :1stofkind, :nthofkind])
    #future_inflation_rate::
    land_space::Entry{T} = Entry(T, "acres", "Plant site space required in acres"; default=1000.0)
    building_volume::Entry{T} = Entry(T, "m^3", "Volume of the tokmak building"; default=140.0e3)
    interest_rate::Entry{T} = Entry(T, "-", "Annual interest rate fraction of direct capital cost"; default=0.05)
    fixed_charge_rate::Entry{T} = Entry(T, "-", "Constant dollar fixed charge rate"; default=0.078)
    indirect_cost_rate::Entry{T} = Entry(T, "-", "Indirect cost associated with construction, equipment, services, energineering construction management and owners cost"; default=0.4)
    lifetime::Entry{Int} = Entry(Int, "years", "lifetime of the plant"; default=40)
    availability::Entry{T} = Entry(T, "-", "availability fraction of the plant"; default=0.803)
    escalation_fraction::Entry{T} = Entry(T, "-", "yearly escalation fraction based on risk assessment"; default=0.05)
    blanket_lifetime::Entry{T} = Entry(T, "years", "lifetime of the blanket"; default=6.8)
    initial_cost_blanket::Entry{T} = Entry(T, "millions of dollars", "cost of initial blanket"; default=200)
    initial_cost_divertor::Entry{T} = Entry(T, "millions of dollars", "cost of initial divertor"; default=8)
    divertor_fluence_lifetime::Entry{T} = Entry(T, "MW/yr*m^2", "divertor fluence lifetime"; default=10)
    blanket_fluence_lifetime::Entry{T} = Entry(T, "MW/yr*m^2", "blanket fluence lifetime"; default=15)
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

    empty!(cst)

    if par.model == :ARIES
        costing_ARIES(dd, par)
    elseif par.model == :Sheffield
        costing_Sheffield(dd, par)
    elseif par.model == :FUSE
        costing_FUSE(dd, par)
    end

    return actor
end

function _finalize(actor::ActorCosting)
    # sort system/subsystem by their costs
    sort!(actor.dd.costing.cost_direct_capital.system, by=x -> x.cost, rev=true)
    for sys in actor.dd.costing.cost_direct_capital.system
        sort!(sys.subsystem, by=x -> x.cost, rev=true)
    end
end

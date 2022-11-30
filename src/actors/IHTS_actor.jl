#= =================== =#
#  ActorIHTS (INTERMEDIATE HEAT TRANSFER SYSTEM) #
#= =================== =#
# ACTOR FOR THE INTERMEDIATE HEAT TRANSFER SYSTEM

mutable struct ActorIHTS <: FacilityAbstractActor
    dd::IMAS.dd
    par::ParametersActor
    tmax
    div_util
    blk_util
end


function ParametersActor(::Type{Val{:ActorIHTS}})
    par = ParametersActor(nothing)
    par.tmax                = Entry(Real, "", "Maximum Outlet Temperature"; default=550)
    par.div_util            = Entry(Real, "", "Utilization Factor of Divertor Circuit"; default=0.85)
    par.blk_util            = Entry(Real, "", "Utilization Factor of Blanket Circuit"; default=0.85)
    return par
end

"""
    ActorIHTS(dd::IMAS.dd, act::ParametersAllActors; kw...)

!!! note 
    Stores data in `dd.balance_of_plant`
"""
function ActorIHTS(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorIHTS(kw...)
    actor = ActorIHTS(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorIHTS(dd::IMAS.dd, par::ParametersActor; kw...)
    logging_actor_init(ActorIHTS)
    par = par(kw...)
    return ActorIHTS(dd, par, par.tmax,par.div_util,par.blk_util)
end


function _step(actor::ActorIHTS)
    dd = actor.dd
    bop = dd.balance_of_plant
    bop.time = dd.core_profiles.time
    # ======= #
    # THERMAL #
    # ======= #
    bop_IHTS = bop.IHTS
    blk_pow_in = [sum([bmod.time_slice[time].power_thermal_extracted for bmod in dd.blanket.module]) for time in bop.time]
    div_pow_in= sum([IMAS.get_time_array(div.power_thermal_extracted, :data, bop.time, :constant) for div in dd.divertors.divertor])
    
    bop_IHTS.divertor_heat_power    = div_pow_in.*actor.div_util
    bop_IHTS.blanket_heat_power     = blk_pow_in.*actor.blk_util
    return actor
end
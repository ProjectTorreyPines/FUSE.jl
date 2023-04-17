#= ================= =#
#  ActorHeatTransfer
#= ================= =#
# ACTOR FOR THE INTERMEDIATE HEAT TRANSFER SYSTEM

const coolant_fluid = [:He, :PbLi]

Base.@kwdef mutable struct FUSEparameters__ActorHeatTransfer{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(Nothing)
    _name::Symbol = :not_set
    #  BREEDER INFO
    breeder_pressure::Entry{T} = Entry(T, "Pa", "Pressure across pump in breeder fluid circuit"; default=1e6)
    breeder_ΔP::Entry{T} = Entry(T, "Pa", "Pressure drop during cooling and heat exchanger"; default=0.25 * 10^6)
    breeder_low_temp::Entry{T} = Entry(T, "K", "Minimum breeder fluid temperature after primary heat exchange"; default=700 + 273.15)
    breeder_hi_temp::Entry{T} = Entry(T, "K", "Maximum breeder fluid temperature at Wall Outlet"; default=1100 + 273.15)
    breeder_η_pump::Entry{T} = Entry(T, "-", "Breeder pump effeciency"; default=0.7)
    #  WALL COOLING INFO
    coolant_pressure::Entry{T} = Entry(T, "Pa", "Pressure across pump in coolant fluid ciruit"; default=10e6)
    coolant_ΔP::Entry{T} = Entry(T, "Pa", "Pressure drop during cooling and heat exchanger"; default=0.3 * 10^6)
    divertor_η_pump::Entry{T} = Entry(T, "-", "Divertor pump effeciency"; default=0.89)
    divertor_max_temp::Entry{T} = Entry(T, "K", "Divertor maximum coolant outlet temperature"; default=650 + 273.15)
    blanket_η_pump::Entry{T} = Entry(T, "-", "Wall pump effeciency"; default=0.89)
    blanket_max_temp::Entry{T} = Entry(T, "K", "Wall maximum coolant outlet temperature"; default=450 + 273.15)
    # ASSUMED 
    breeder_HX_ϵ::Entry{T} = Entry(T, "-", "Effectiveness of the breeder - cycle heat exchanger"; default=0.9)
    divertor_HX_ϵ::Entry{T} = Entry(T, "-", "Effectiveness of the divertor - cycle heat exchanger"; default=0.9)
    blanket_HX_ϵ::Entry{T} = Entry(T, "-", "Effectiveness of the wall - cycle heat exchanger"; default=0.9)
    breeder_fluid::Switch{Symbol} = Switch(Symbol, coolant_fluid, "-", "Breeder coolant fluid"; default=:PbLi)
    blanket_coolant::Switch{Symbol} = Switch(Symbol, coolant_fluid, "-", "Breeder coolant fluid"; default=:He)
    divertor_coolant::Switch{Symbol} = Switch(Symbol, coolant_fluid, "-", "Breeder coolant fluid"; default=:He)
end

mutable struct ActorHeatTransfer <: FacilityAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorHeatTransfer
end

"""
    ActorHeatTransfer(dd::IMAS.dd, act::ParametersAllActors; kw...)

!!! note 
    Stores data in `dd.balance_of_plant`
"""
function ActorHeatTransfer(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorHeatTransfer
    actor = ActorHeatTransfer(dd, par; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorHeatTransfer(dd::IMAS.dd, par::FUSEparameters__ActorHeatTransfer, act::ParametersAllActors; kw...)
    logging_actor_init(ActorHeatTransfer)
    par = par(kw...)
    return ActorHeatTransfer(dd, par)
end

function _step(actor::ActorHeatTransfer)
    dd = actor.dd
    par = actor.par
    bop = dd.balance_of_plant

    bop.heat_transfer.divertor.working_fluid = string(par.divertor_coolant)
    bop.heat_transfer.wall.working_fluid = string(par.blanket_coolant)
    bop.heat_transfer.breeder.working_fluid = string(par.breeder_fluid)

    bop.time = dd.core_profiles.time

    # ======= #
    # THERMAL #
    # ======= #
    bop_IHTS = bop.heat_transfer

    breeder_heat_load = 0.0
    if !isempty(dd.blanket.module)
        breeder_heat_load = sum(bmod.time_slice[].power_thermal_extracted for bmod in dd.blanket.module)
    end
    divertor_heat_load = 0.0
    if !isempty(dd.divertors.divertor)
        divertor_heat_load = sum((@ddtime(div.power_incident.data)) for div in dd.divertors.divertor)
    end
    blanket_heat_load = abs.(IMAS.radiation_losses(dd.core_sources))

    # @show breeder_heat_load
    # @show divertor_heat_load
    # @show blanket_heat_load

    cp_low, rho_low = pbLi_props(par.breeder_low_temp)
    cp_hi, rho_hi = pbLi_props(par.breeder_hi_temp)
    cp_ave = (cp_hi + cp_low) / 2.0
    rho_ave = (rho_hi + rho_low) / 2.0
    v_ave = 1.0 / rho_ave

    cp_he = 5.1926e3

    cv_cycle = 3.1156e3
    k_cycle = cp_he / cv_cycle
    kcoeff = (k_cycle - 1.0) / k_cycle

    # inital values as guess
    ΔT_cycle_breeder = par.breeder_hi_temp - par.divertor_max_temp
    cycle_flow = breeder_heat_load / (ΔT_cycle_breeder * cp_he)
    @ddtime(bop.thermal_cycle.flow_rate = cycle_flow)
    @ddtime(bop_IHTS.wall.flow_rate = cycle_flow)
    @ddtime(bop_IHTS.divertor.flow_rate = cycle_flow)

    # INTIALIZING DD WITH HEAT LAODS
    @ddtime(bop_IHTS.divertor.heat_load = divertor_heat_load)
    @ddtime(bop_IHTS.wall.heat_load = blanket_heat_load)
    @ddtime(bop_IHTS.breeder.heat_load = breeder_heat_load)

    #basic params
    ΔT_blanket = blanket_heat_load / (cycle_flow * cp_he)
    ΔT_divertor = divertor_heat_load / (cycle_flow * cp_he)

    rp_div = par.coolant_pressure / (par.coolant_pressure - par.coolant_ΔP)
    rp_blk = rp_div

    @ddtime(bop_IHTS.wall.outlet_temperature = par.blanket_max_temp)
    @ddtime(bop_IHTS.divertor.outlet_temperature = par.divertor_max_temp)

    blk_after_comp = par.blanket_max_temp - ΔT_blanket
    div_after_comp = par.divertor_max_temp - ΔT_divertor

    @ddtime(bop_IHTS.wall.inlet_temperature = blk_after_comp)
    @ddtime(bop_IHTS.divertor.inlet_temperature = div_after_comp)

    blk_pre_comp = blk_after_comp / (1.0 + (rp_blk^kcoeff - 1.0) / par.blanket_η_pump)
    div_pre_comp = div_after_comp / (1.0 + (rp_div^kcoeff - 1.0) / par.divertor_η_pump)

    @ddtime(bop_IHTS.wall.HX_outlet_temperature = blk_pre_comp)
    @ddtime(bop_IHTS.divertor.HX_outlet_temperature = div_pre_comp)

    @ddtime(bop_IHTS.wall.circulator_power = gas_circulator(rp_blk, blk_pre_comp, par.blanket_η_pump, cycle_flow))
    @ddtime(bop_IHTS.divertor.circulator_power = gas_circulator(rp_div, div_pre_comp, par.divertor_η_pump, cycle_flow))

    w_isentropic = v_ave * par.breeder_ΔP
    w_actual = w_isentropic / par.breeder_η_pump

    Tin_breeder = par.breeder_low_temp + (w_actual - w_isentropic) / cp_ave   #temperature after pump
    ΔT_breeder = par.breeder_hi_temp - Tin_breeder
    breeder_mflow = breeder_heat_load / (cp_ave * ΔT_breeder)

    # FULLY DEFINING BREEDER
    @ddtime(bop_IHTS.breeder.HX_outlet_temperature = par.breeder_low_temp)
    @ddtime(bop_IHTS.breeder.inlet_temperature = Tin_breeder)
    @ddtime(bop_IHTS.breeder.outlet_temperature = par.breeder_hi_temp)
    @ddtime(bop_IHTS.breeder.flow_rate = breeder_mflow)
    @ddtime(bop_IHTS.breeder.circulator_power = breeder_mflow * v_ave * par.breeder_ΔP / par.breeder_η_pump)

    return actor
end

function gas_circulator(rp, Tin, effC, mflow, nstages=1)
    cp = 5.1926e3
    cv = 3.1156e3
    a1c = (effC * (rp^(1.0 / nstages))^(cv / cp) - (rp^(1.0 / nstages))^(cv / cp) + rp^(1.0 / nstages)) / (effC * (rp^(1.0 / nstages))^(cv / cp))
    work_in = cp * (a1c - 1.0) * Tin * mflow
    return work_in
end
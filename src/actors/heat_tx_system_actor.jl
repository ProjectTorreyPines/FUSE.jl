#= ================= =#
#  ActorHeatTxSystem
#= ================= =#
# ACTOR FOR THE INTERMEDIATE HEAT TRANSFER SYSTEM

mutable struct ActorHeatTxSystem<: FacilityAbstractActor
    dd::IMAS.dd
    par::ParametersActor
    radiation_factor::Float64
end


function ParametersActor(::Type{Val{:ActorHeatTxSystem}})
    par = ParametersActor(nothing)
    #  BREEDER INFO
    par.breeder_pressure     = Entry(Real, "", "Pressure across pump in breeder fluid circuit"; default = 1e6)
    par.breeder_ΔP          = Entry(Real,"",    "Pressure drop during cooling and heat exchanger"; default = 0.25*10^6)
    par.breeder_low_temp    = Entry(Real, "", "Minimum breeder fluid temperature after primary heat exchange"; default=700+273.15)
    par.breeder_hi_temp     = Entry(Real, "", "Maximum breeder fluid temperature at Blanket Outlet"; default=1100+273.15)
    par.breeder_η_pump      = Entry(Real, "", "Breeder pump effeciency"; default=0.7)
    #  BLANKET COOLING INFO
    par.coolant_pressure        = Entry(Real, "",   "Pressure across pump in coolant fluid ciruit"; default = 10e6)
    par.coolant_ΔP              = Entry(Real,"",    "Pressure drop during cooling and heat exchanger"; default = 0.3*10^6)
    
    par.divertor_η_pump         = Entry(Real, "",  "Divertor pump effeciency"; default=0.89)
    par.divertor_max_temp       = Entry(Real, "",   "Maximum Coolant Outlet Temperature"; default=650+273.15)
    
    par.blanket_η_pump          = Entry(Real, "",   "Blanket pump effeciency"; default=0.89)
    par.blanket_max_temp        = Entry(Real, "",   "Maximum Coolant Outlet Temperature"; default=450+273.15)


    # ASSUMED 
    par.radiation_factor    = Entry(Real, "",   "Assumed factor of multiplication for the power absorbed by the blanket";default=2.5)
    par.breeder_HX_ϵ        = Entry(Real, "",   "effectiveness of the breeder - cycle heat exchanger";default=0.9)
    par.divertor_HX_ϵ       = Entry(Real, "",   "effectiveness of the divertor - cycle heat exchanger";default=0.9)
    par.blanket_HX_ϵ        = Entry(Real, "",   "effectiveness of the blanket - cycle heat exchanger";default=0.9)
    par.breeder_fluid       = Entry(String,"",  "Divertor coolant fluid";default = "pbli")
    par.blanket_coolant     = Entry(String,"",  "Divertor coolant fluid";default = "helium")
    par.divertor_coolant    = Entry(String,"",  "Divertor coolant fluid";default = "helium")
    return par
end

"""
    ActorHeatTxSystem(dd::IMAS.dd, act::ParametersAllActors; kw...)

!!! note 
    Stores data in `dd.balance_of_plant`
"""

function ActorHeatTxSystem(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorHeatTxSystem(kw...)
    bop = dd.balance_of_plant;
    bop.heat_tx_system.divertor.working_fluid  = par.divertor_coolant
    bop.heat_tx_system.blanket.working_fluid   = par.blanket_coolant
    bop.heat_tx_system.breeder.working_fluid   = par.breeder_fluid
    actor = ActorHeatTxSystem(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorHeatTxSystem(dd::IMAS.dd, par::ParametersActor; kw...)
    logging_actor_init(ActorHeatTxSystem)
    par = par(kw...)
    return ActorHeatTxSystem(dd, par, par.radiation_factor)
end


function _step(actor::ActorHeatTxSystem)
    dd = actor.dd
    par = actor.par
    bop = dd.balance_of_plant
    bop.time = dd.core_profiles.time
    # ======= #
    # THERMAL #
    # ======= #
    bop_IHTS = bop.heat_tx_system
    
    # breeder_heat_load   = (sum([bmod.time_slice[time].power_thermal_extracted for bmod in dd.blanket.module]) for time in bop.time)
    breeder_heat_load   = abs.(sum([bmod.time_slice[].power_thermal_extracted for bmod in dd.blanket.module]))
    divertor_heat_load  = abs.(sum([(@ddtime(div.power_incident.data)) for div in dd.divertors.divertor]))
    blanket_heat_load   = abs.(IMAS.radiation_losses(dd.core_sources)*actor.radiation_factor);
    
    # @show divertor_heat_load
    # @show blanket_heat_load
    # @show breeder_heat_load

    cp_low,rho_low  = pbLi_props(par.breeder_low_temp);
    cp_hi,rho_hi    = pbLi_props(par.breeder_hi_temp);
    cp_ave = (cp_hi+cp_low)/2;
    rho_ave = (rho_hi+rho_low)/2;
    v_ave = 1/rho_ave;

    cp_he = 5.1926e3;

    cv_cycle    = 3.1156e3;    
    k_cycle = cp_he/cv_cycle;
    kcoeff      = (k_cycle-1)/k_cycle;

    # inital values as guess
    ΔT_cycle_breeder = par.breeder_hi_temp - par.divertor_max_temp;

    @ddtime(bop.thermal_cycle.flow_rate = breeder_heat_load/(ΔT_cycle_breeder*cp_he))
    
    cycle_flow = @ddtime(bop.thermal_cycle.flow_rate) ;

    @ddtime(bop_IHTS.blanket.flow_rate = cycle_flow);
    @ddtime(bop_IHTS.divertor.flow_rate = cycle_flow);

        # INTIALIZING DD WITH HEAT LAODS
    @ddtime(bop_IHTS.divertor.heat_load    = divertor_heat_load)
    @ddtime(bop_IHTS.blanket.heat_load     = blanket_heat_load)
    @ddtime(bop_IHTS.breeder.heat_load     = breeder_heat_load)

    #basic params
    ΔT_blanket  = blanket_heat_load/(cycle_flow*cp_he)
    ΔT_divertor = divertor_heat_load/(cycle_flow*cp_he)
    
    rp_div = par.coolant_pressure/(par.coolant_pressure-par.coolant_ΔP);
    rp_blk = rp_div;

    @ddtime(bop_IHTS.blanket.outlet_temperature = par.blanket_max_temp)
    @ddtime(bop_IHTS.divertor.outlet_temperature = par.divertor_max_temp)

    blk_after_comp =  par.blanket_max_temp-ΔT_blanket
    div_after_comp = par.divertor_max_temp-ΔT_divertor
    
    @ddtime(bop_IHTS.blanket.inlet_temperature  = blk_after_comp)
    @ddtime(bop_IHTS.divertor.inlet_temperature = div_after_comp)

    blk_pre_comp = blk_after_comp/(1+(rp_blk^kcoeff-1)/par.blanket_η_pump);
    div_pre_comp = div_after_comp/(1+(rp_div^kcoeff-1)/par.divertor_η_pump);

    @ddtime(bop_IHTS.blanket.HX_outlet_temperature  =   blk_pre_comp)
    @ddtime(bop_IHTS.divertor.HX_outlet_temperature = div_pre_comp)

    @ddtime(bop_IHTS.blanket.circulator_power  =    gasCirculator(rp_blk,blk_pre_comp,par.blanket_η_pump,cycle_flow))
    @ddtime(bop_IHTS.divertor.circulator_power  =   gasCirculator(rp_div,div_pre_comp,par.divertor_η_pump,cycle_flow))

    w_isentropic    = v_ave*par.breeder_ΔP
    w_actual        = w_isentropic/par.breeder_η_pump

    Tin_breeder =   par.breeder_low_temp+(w_actual-w_isentropic)/cp_ave;   #temperature after pump
    ΔT_breeder  =   par.breeder_hi_temp-Tin_breeder;
    breeder_mflow =  breeder_heat_load/(cp_ave*ΔT_breeder);

        # FULLY DEFINING BREEDER
    @ddtime(bop_IHTS.breeder.HX_outlet_temperature  = par.breeder_low_temp);
    @ddtime(bop_IHTS.breeder.inlet_temperature  = Tin_breeder);
    @ddtime(bop_IHTS.breeder.outlet_temperature = par.breeder_hi_temp)
    @ddtime(bop_IHTS.breeder.flow_rate          = breeder_mflow);
    @ddtime(bop_IHTS.breeder.circulator_power   = breeder_mflow*v_ave*par.breeder_ΔP/par.breeder_η_pump);

    return actor
end

function pbLi_props(Temperature)
#temperature input in celcius
specific_heat   = (0.195 - 9.116*10^(-6).*(Temperature)).*1000   #J/kgK
density         =  10520.35-1.19051.*(Temperature); 
return [specific_heat,density]
end

function gasCirculator(rp,Tin,effC,mflow)
    nstages = 1;

    cp = 5.1926e3;
    cv = 3.1156e3;
    a1c = (effC*(rp^(1/nstages))^(cv/cp) - (rp^(1/nstages))^(cv/cp) + rp^(1/nstages))/(effC*(rp^(1/nstages))^(cv/cp));

    work_in   = cp*(a1c-1)*Tin*mflow;
    return work_in
end
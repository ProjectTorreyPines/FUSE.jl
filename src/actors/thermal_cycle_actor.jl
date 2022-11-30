#= =================== =#
#  Thermal Cycle actor  #
#= =================== =#
using LinearAlgebra
using Plots

# CONSTRUCTOR - SHOULD Hold 
mutable struct ActorThermalCycle <: FacilityAbstractActor
    dd::IMAS.dd
    par::ParametersActor
    cyc_name::String
end

function ParametersActor(::Type{Val{:ActorThermalCycle}})
    par         = ParametersActor(nothing)
    par.model   = Entry(String, "", "Power Cycle model"; default="brayton")
    par.rp      = Entry(Real, "", "Overall Compression Ratio"; default=3.0)
    par.Pmax    = Entry(Real, "", "Max System Pressure (kPa)"; default=8000)
    par.Tmax    = Entry(Real,"","Max Cycle Temperature K";default = 600.0+273.15)
    par.Tmin    = Entry(Real,"","Min Cycle Temperature K";default = 35.0+273.15)
    par.Nt      = Entry(Int, "", "Number of Turbine Stages"; default=2)
    par.Nc      = Entry(Int, "", "Number of Compression Stages"; default=2)
    return par
end
"""
    ActorThermalCYcle(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates the thermal cycle based on the heat loads of the divertor and the blanket
* `model = :brayton` evaluates the brayton cycle
* `model = :rankine evaluates the rankine cycle **THIS IS NOT ADDED YET, WILL UPDATE SOON
!!! note 
    Stores data in `dd.balance_of_plant`
"""
function ActorThermalCycle(dd::IMAS.dd, act::ParametersAllActors; kw...)    #constructor 2
    par = act.ActorThermalCycle(kw...)                                      #CALLING ParametersActor
    actor = ActorThermalCycle(dd, par)                                      #Calling constructor 3, returns actor with all parameters
    step(actor)
    finalize(actor)
    return actor
end

function ActorThermalCycle(dd::IMAS.dd, par::ParametersActor; kw...)            #constructor 3
    logging_actor_init(ActorThermalCycle)
    par = par(kw...)
    return ActorThermalCycle(dd, par,par.model)                                          #CALLING CONTRUCTOR ON LINE 9
end


function _step(actor::ActorThermalCycle)
    dd = actor.dd
    bop = dd.balance_of_plant

    # ======= #
    # THERMAL #
    # ======= #
    bop_thermal = bop.thermal_cycle

    ### BLANKET AND DIVERTOR POWER###
    # blanket_power   =       [sum([bmod.time_slice[time].power_thermal_extracted for bmod in dd.blanket.module]) for time in bop.time]
    # div_power       =       sum([IMAS.get_time_array(div.power_thermal_extracted, :data, bop.time, :constant) for div in dd.divertors.divertor])
    div_power       = abs.(bop.IHTS.divertor_heat_power)
    blanket_power   = abs.(bop.IHTS.blanket_heat_power)

    bop_thermal.total_heat_power = div_power+blanket_power  #power from divertor and blanket
    braytonOut = braytonCycle(actor.par.rp,actor.par.Pmax,actor.par.Tmin,actor.par.Tmax, actor.par.Nt, actor.par.Nc)

    bop_thermal.thermal_effeciency      = braytonOut.η_th.*ones(length(dd.core_profiles.time))
    bop_thermal.mass_flow_rate          = bop_thermal.total_heat_power./braytonOut.η_th
    bop_thermal.cycle_work_output       = bop_thermal.mass_flow_rate .* braytonOut.w_out
    bop_thermal.cycle_work_input        = bop_thermal.mass_flow_rate .* braytonOut.w_in
    bop_thermal.cycle_net_work          = bop_thermal.mass_flow_rate .* braytonOut.w_net

    bop_thermal.heat_waste              = bop_thermal.mass_flow_rate .* braytonOut.q_L
    bop_thermal.IHTS_inlet_temperature  = braytonOut.T_HX .* ones(length(dd.core_profiles.time))

    return actor
end

struct BraytonOutput{T<:Real}
    η_th::T
    w_net::T
    w_out::T
    w_in::T
    r_bw::T
    q_H ::T
    q_L ::T
    T_HX::T
end
function braytonCycle(rp::Real,Pmax::Real,Tmin::Real,Tmax::Real, Nt::Int, Nc::Int; ηt::Real=0.93, ηc::Real=0.89, ϵr::Real=0.95)
    #Brayton evaluates the thermal performance for a specified brayton cycle
    #   
    # ===================================================================
    # INPUTS        UNIT        DESCRIPTION
    #   rp          Scalar,     Total Compression Ratio
    #   Pmax        [kPa],      Specified pressure, maximum
    #   Tmin        [K]         Lowest Cycle Temp
    #   Tmax        [K]         Highest Cycle Temp
    #   Nc,Nt       Scalars     Number of compression (Nc) and Turbine (Nt) stages
    #   ηc          Frac        Fractional compressor isentropic effeciency
    #   ηt          Frac        Fractional Turbine isentropic effeciency
    #   ϵr          Frac        Fractional Regenerator effectiveness

    cp = 5.1926;
    cv = 3.1156;
    
    #SUMMARIZED COEFFECICENTS
    a1t = (ηt*(rp^(1/Nt))^(cv/cp) - ηt*rp^(1/Nt) + rp^(1/Nt))/rp^(1/Nt)
    a1c = (ηc*(rp^(1/Nc))^(cv/cp) - (rp^(1/Nc))^(cv/cp) + rp^(1/Nc))/(ηc*(rp^(1/Nc))^(cv/cp))
    Pmin = Pmax/rp
    
    #STATE VARIABLES
    phi = [Tmin;Tmax]
    θ_L = [Tmin;Pmin]
    θ_H = [Tmax;Pmax]
    
   #MATRICES, DEFINED IN WRITE UP
    A_t = Diagonal([a1t,rp^(1/Nt)])    
    B_t = Diagonal([1/a1t,1])           #turbine
    A_c = Diagonal([a1c,rp^(-1/Nc)])
    B_c = Diagonal([1/a1c,1])#Comp
    C   = Diagonal([cp*(a1c-1),cp*(a1t-1)])
    N   = [Nc 0;(1-Nc) 0;0 Nt;0 (1-Nt)]
    E   = [(1-ϵr) ϵr;ϵr (1-ϵr)]
    
    #OUTPUTS OF THE COMPRESSOR AND TURBINE CIRCUITS
    θ_ci = (A_c*B_c)^(Nc-1)*A_c*θ_L    #AFTER COMP
    θ_ti = (A_t*B_t)^(Nt-1)*A_t*θ_H    #AFTER TURB

    θ_1 = E*[θ_ci[1];θ_ti[1]]  #OUTLET TEMPERATURES OF REGENERATOR
    θ_2 = [Tmax;Tmin]          #FOR USE IN THE FOLLOWING

    #RESULTS
    wc,q_intercool,wt,q_reheat  = N*C*phi;                      #[COMP WORK, INTERCOOL HEAT LOSS, TURBINE WORK, REHEAT HEAD ADDITION]
    q_primary,q_waste           = cp*(θ_2-θ_1);                 #FUSION HEAT ADDITION, HEAT WASTED
    wnet                        = abs(wt)-abs(wc);              #NET WORK OUT 
    q_in                        = abs(q_primary)+abs(q_reheat); #NET HEAT ADDED TO THE SYSTEM
    q_out                       = abs(q_waste)+abs(q_intercool);#NET HEAT REMOVED FROM SYSTEM
    η_thermal                   = wnet/q_in;                    #THERMAL EFFECIENCY
    Primary_inlet_temp          = θ_1[1];                       #INLET HELIUM TEMPERATURE AT REACTOR HEAT EXCHANGERS

    # ALL  return result = (inputs = [rp,Pmax,Tmin,Tmax,Nt,Nc,ηt,ηc,ϵr], netWork = wnet, turbWork = wt, compWork = wc, heatIn = q_in, heatWaste = q_waste+q_intercool, η_th = η_thermal)
    # SOME return result = (netWork = wnet, turbWork = wt, compWork = wc, heatIn = q_in, heatWaste = q_waste+q_intercool, η_th = η_thermal)
    # SOME return result = (wc = wc,wt=wt,q_rh = q_reheat,q_ic = q_intercool,q_fus = q_primary,q_rej = q_waste,T_HX = Primary_inlet_temp)
    # ONLY NECECCARY
    return BraytonOutput(η_thermal, abs(wnet),abs(wc),abs(wt), abs(wc/wt), abs(q_in) ,abs(q_out),Primary_inlet_temp)
end
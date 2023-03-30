#= =================== =#
#  Thermal Cycle actor  #
#   Cases
#       SIMPLE Brayton
#       COMPLEX BRAYTON
#       SIMPLE rankine
#       COMBINED BRAYTON RANKINE
#= =================== =#

Base.@kwdef mutable struct FUSEparameters__ActorThermalCycle{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(Nothing)
    _name::Symbol = :not_set
    power_cycle_type::Switch{Symbol} = Switch(Symbol, [:brayton_only, :rankine_only, :complex_brayton], "-", "Power cycle configuration"; default=:complex_brayton)
    rp::Entry{T} = Entry(T, "-", "Overall compression ratio"; default=3.0)
    Pmax::Entry{T} = Entry(T, "-", "Max system pressure (MPa)"; default=8e6)
    Tmax::Entry{T} = Entry(T, "-", "Max cycle temperature K"; default=950.0 + 273.15)
    Tmin::Entry{T} = Entry(T, "-", "Min cycle temperature K"; default=35.0 + 273.15)
    Nt::Entry{Int} = Entry(Int, "-", "Number of turbine stages"; default=1)
    Nc::Entry{Int} = Entry(Int, "-", "Number of compression stages"; default=3)
    regen::Entry{T} = Entry(T, "-", "Regeneration fraction")
    do_plot::Entry{Bool} = Entry(Bool, "-", "plot"; default=false)
end

mutable struct ActorThermalCycle <: FacilityAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorThermalCycle
    act::ParametersAllActors
    function ActorThermalCycle(dd::IMAS.dd, par::FUSEparameters__ActorThermalCycle, act::ParametersAllActors; kw...)
        logging_actor_init(ActorThermalCycle)
        par = par(kw...)
        return new(dd, par, act)
    end
end

"""
    ActorThermalCycle(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates the thermal cycle based on the heat loads of the divertor and the wall
* `model = :brayton` evaluates the brayton cycle
* `model = :rankine evaluates the rankine cycle
!!! note 
    Stores data in `dd.balance_of_plant`
"""
function ActorThermalCycle(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorThermalCycle(kw...)
    actor = ActorThermalCycle(dd, par, act)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorThermalCycle)
    dd = actor.dd
    par = actor.par
    bop = dd.balance_of_plant
    ihts = bop.heat_transfer
    wall = ihts.wall
    divertor = ihts.divertor

    bop = dd.balance_of_plant
    bop.power_cycle_type = string(par.power_cycle_type)

    bop_thermal = bop.thermal_cycle
    ihts_par = actor.act.ActorHeatTransfer

    blanket_power = @ddtime(bop.heat_transfer.wall.heat_load)
    breeder_power = @ddtime(bop.heat_transfer.breeder.heat_load)
    divertor_power = @ddtime(bop.heat_transfer.divertor.heat_load)

    if ismissing(par, :regen)
        ϵr = 0.0
        if bop.power_cycle_type ∈ ["brayton_only"]
            ϵr = 0.9
        end
    else
        ϵr = par.regen
    end
    @assert (ϵr == 0 || (ϵr >= 0.0 && ϵr <= 1.0 && bop.power_cycle_type ∈ ["brayton_only"])) "Regeneration between 0.0 and 1.0 is only possible for `:brayton_only` cycles"

    mflow_cycle = @ddtime(bop.thermal_cycle.flow_rate)

    if bop.power_cycle_type == "rankine_only"
        return actor

    elseif bop.power_cycle_type == "complex_brayton"
        cp_cycle = 5.1926e3
        cv_cycle = 3.1156e3
        k_cycle = cp_cycle / cv_cycle
        kcoeff = (k_cycle - 1) / k_cycle
        cp_blk = cp_cycle
        kcoeff_blk = kcoeff
        cp_div = cp_cycle
        kcoeff_div = kcoeff

        rp_blanket = ihts_par.coolant_pressure / (ihts_par.coolant_pressure - ihts_par.coolant_ΔP)
        rp_divertor = rp_blanket

        Niter = 100
        cp_he = 5.1926e3
        mratio = 1.0

        pc = compressor_stages(par)

        @ddtime(bop_thermal.input_work = mflow_cycle * pc.work)

        hxA = ihts_heat_exchanger(ihts_par.blanket_η_pump, ihts_par.blanket_HX_ϵ, Niter, pc.T_out, @ddtime(bop.thermal_cycle.flow_rate), cp_he, mratio, blanket_power, rp_blanket, cp_blk, kcoeff_blk)

        @ddtime(wall.HX_outlet_temperature = hxA.Tmin)
        @ddtime(wall.inlet_temperature = hxA.Tin)
        @ddtime(wall.outlet_temperature = hxA.Tmax)
        @ddtime(wall.circulator_power = hxA.pump_work)
        @ddtime(wall.heat_delivered = hxA.HX_q)
        @ddtime(wall.heat_waste = abs(blanket_power + hxA.pump_work - hxA.HX_q))
        hxB = ihts_heat_exchanger(ihts_par.divertor_η_pump, ihts_par.divertor_HX_ϵ, Niter, hxA.Tout_cycle, @ddtime(bop.thermal_cycle.flow_rate), cp_he, mratio, divertor_power, rp_divertor, cp_div, kcoeff_div)

        @ddtime(divertor.HX_outlet_temperature = hxB.Tmin)
        @ddtime(divertor.inlet_temperature = hxB.Tin)
        @ddtime(divertor.outlet_temperature = hxB.Tmax)
        @ddtime(divertor.circulator_power = hxB.pump_work)
        @ddtime(divertor.heat_delivered = hxB.HX_q)
        @ddtime(divertor.heat_waste = divertor_power + hxB.pump_work - hxB.HX_q)

        Tb_max = @ddtime(bop.heat_transfer.breeder.outlet_temperature)
        Tb_min = @ddtime(bop.heat_transfer.breeder.HX_outlet_temperature)
        Tb_blanket_in = @ddtime(bop.heat_transfer.breeder.inlet_temperature)
        Tbreeder_ave = (Tb_max + Tb_min) / 2
        cp_pbli, rho_pbli = pbLi_props(Tbreeder_ave)
        breeder_pump_ΔT = Tb_blanket_in - Tb_min

        if abs(par.Tmax - Tb_max) > 35.0
            @debug "ActorThermalCycle: adjusting preset max cycle temperature"
            actor.par.Tmax = Tb_max - 50.0
            par = actor.par
        end

        mflow_breeder = @ddtime(bop.heat_transfer.breeder.flow_rate)
        bpump = @ddtime(bop.heat_transfer.breeder.circulator_power)

        mcp_breeder = cp_pbli * mflow_breeder
        mcp_cyc = mflow_cycle * cp_he
        OG_tbmin = Tb_min
        Tb_after_hx = Tb_min
        post_turbine_guess = turbine_stages(par)

        qq = 0.0
        Tco = 0.0
        Tcyc_out = 0.0
        #iterating

        Tci_regen_cycle = hxB.Tout_cycle
        ηt = 0.93
        ϵ_breeder = ihts_par.breeder_HX_ϵ
        ϵ_regen = 0.9
        turb_coeff = 1 - ηt + (1 / par.rp)^(kcoeff) * ηt
        Cmin = minimum([mcp_breeder, mcp_cyc])
        C_coeff = ϵ_breeder * (Cmin / mcp_cyc)

        A = [1.0 0.0 -ϵ_regen; C_coeff 1.0 0.0; 0.0 -turb_coeff 1.0]
        b = [(Tci_regen_cycle * (1.0 - ϵ_regen)); C_coeff * Tb_max; 0.0]

        x = A \ b
        par.Tmax = x[2]

        for iter in 1:50
            post_turbine_guess = turbine_stages(par)
            Tco, Tho = regenHX(hxB.Tout_cycle, post_turbine_guess.T_out, 0.95)

            qq, Tbreed_out, Tcyc_out = hxeff(Tb_max, mcp_breeder, Tco, mcp_cyc, ihts_par.breeder_HX_ϵ)

            par.Tmax = Tcyc_out
            Tb_after_hx = Tbreed_out

            if Tb_after_hx >= (OG_tbmin - 5.0)
                break
            end

            Tb_min = Tb_after_hx + 5.0
            Tb_blanket_in = Tb_min + breeder_pump_ΔT
            ΔT_breeder_blanket = Tb_max - Tb_min

            #something is wrong here, i need to fix it
            oldMflow = mflow_breeder
            mflow_breeder = breeder_power / (ΔT_breeder_blanket * cp_pbli)
            bpump = bpump * (mflow_breeder / oldMflow)
            mcp_breeder = mflow_breeder * cp_pbli

            # ΔT_breeder_blanket  = breeder_power/(mflow_breeder*cp_pbli)
            Tb_max = Tb_blanket_in + breeder_power / (mflow_breeder * cp_pbli)
        end
        if Tb_after_hx < (OG_tbmin - 5.0)
            @warn "ActorThermalCycle has not converged"
        end

        pt = turbine_stages(par)
        @ddtime(bop.thermal_cycle.turbine_work = mflow_cycle * pt.work)
        @ddtime(bop.heat_transfer.breeder.HX_outlet_temperature = Tb_after_hx)
        @ddtime(bop.heat_transfer.breeder.excess_temperature = Tb_after_hx - Tb_min)
        @ddtime(bop.heat_transfer.breeder.inlet_temperature = Tb_blanket_in)
        @ddtime(bop.heat_transfer.breeder.outlet_temperature = Tb_max)
        @ddtime(bop.heat_transfer.breeder.heat_delivered = qq)
        @ddtime(bop.heat_transfer.breeder.heat_waste = breeder_power + bpump - qq)
        @ddtime(bop.thermal_cycle.net_work = mflow_cycle * (pt.work - pc.work))
        @ddtime(bop_thermal.thermal_effeciency = mflow_cycle * (pt.work - pc.work) / (hxB.HX_q + hxA.HX_q + qq))
        return actor

    elseif bop.power_cycle_type == "brayton_only"
        braytonT = braytonCycle(actor.par.rp, actor.par.Pmax, actor.par.Tmin, actor.par.Tmax, actor.par.Nt, actor.par.Nc; ϵr)
        par.Tmax = evalBrayton(braytonT, ihts_par, dd)
        cp_cycle = 5.1926e3
        braytonOut = braytonCycle(actor.par.rp, actor.par.Pmax, actor.par.Tmin, actor.par.Tmax, actor.par.Nt, actor.par.Nc; ϵr)
        @ddtime(bop_thermal.turbine_work = mflow_cycle .* braytonOut.w_out)
        @ddtime(bop_thermal.input_work = mflow_cycle .* braytonOut.w_in)
        @ddtime(bop_thermal.net_work = mflow_cycle .* (braytonOut.w_out - braytonOut.w_in))
        @ddtime(bop_thermal.thermal_effeciency = braytonOut.η_th)
        return actor
    else
        error("Power cycle type `$(bop.power_cycle_type)` not recognized")
    end
end

struct BraytonOutput{T<:Real}
    η_th::T
    w_net::T
    w_in::T
    w_out::T
    r_bw::T
    q_H::T
    q_L::T
    T_HX::T
end

"""
    braytonCycle(rp::T, Pmax::T, Tmin::T, Tmax::T, Nt::Int, Nc::Int; ηt::T=0.93, ηc::T=0.89, ϵr::T=0) where {T<:Real}

Brayton evaluates the thermal performance for a specified brayton cycle
      
    ===================================================================
    INPUTS        UNIT        DESCRIPTION
      rp          Scalar,     Total Compression Ratio
      Pmax        [kPa],      Specified pressure, maximum
      Tmin        [K]         Lowest Cycle Temp
      Tmax        [K]         Highest Cycle Temp
      Nc,Nt       Scalars     Number of compression (Nc) and Turbine (Nt) stages
      ηc          Frac        Fractional compressor isentropic effeciency
      ηt          Frac        Fractional Turbine isentropic effeciency
      ϵr          Frac        Fractional Regenerator effectiveness
"""
function braytonCycle(rp::T, Pmax::T, Tmin::T, Tmax::T, Nt::Int, Nc::Int; ηt::T=0.93, ηc::T=0.89, ϵr::T=0) where {T<:Real}
    cp = 5.1926e3
    cv = 3.1156e3

    #SUMMARIZED COEFFECICENTS
    a1t = (ηt * (rp^(1.0 / Nt))^(cv / cp) - ηt * rp^(1.0 / Nt) + rp^(1.0 / Nt)) / rp^(1.0 / Nt)
    a1c = (ηc * (rp^(1.0 / Nc))^(cv / cp) - (rp^(1.0 / Nc))^(cv / cp) + rp^(1.0 / Nc)) / (ηc * (rp^(1.0 / Nc))^(cv / cp))
    Pmin = Pmax / rp

    #STATE VARIABLES
    phi = [Tmin; Tmax]
    θ_L = [Tmin; Pmin]
    θ_H = [Tmax; Pmax]

    #MATRICES, DEFINED IN WRITE UP
    A_t = Diagonal([a1t, rp^(-1.0 / Nt)])
    B_t = Diagonal([1 / a1t, 1.0])           #turbine
    A_c = Diagonal([a1c, rp^(1.0 / Nc)])
    B_c = Diagonal([1 / a1c, 1.0])#Comp
    C = Diagonal([cp * (a1c - 1.0), cp * (a1t - 1.0)])
    N = [Nc 0; (1-Nc) 0; 0 Nt; 0 (1-Nt)]
    E = [(1.0-ϵr) ϵr; ϵr (1.0-ϵr)]

    #OUTPUTS OF THE COMPRESSOR AND TURBINE CIRCUITS
    θ_ci = (A_c * B_c)^(Nc - 1) * A_c * θ_L    #AFTER COMP
    θ_ti = (A_t * B_t)^(Nt - 1) * A_t * θ_H    #AFTER TURB

    θ_1 = E * [θ_ci[1]; θ_ti[1]]  #OUTLET TEMPERATURES OF REGENERATOR
    θ_2 = [Tmax; Tmin]            #FOR USE IN THE FOLLOWING

    #RESULTS
    wc, q_intercool, wt, q_reheat = N * C * phi                      #[COMP WORK, INTERCOOL HEAT LOSS, TURBINE WORK, REHEAT HEAD ADDITION]
    q_primary, q_waste = cp * (θ_2 - θ_1)                 #FUSION HEAT ADDITION, HEAT WASTED
    wnet = abs(wt) - abs(wc)              #NET WORK OUT 
    q_in = abs(q_primary) + abs(q_reheat) #NET HEAT ADDED TO THE SYSTEM
    q_out = abs(q_waste) + abs(q_intercool)#NET HEAT REMOVED FROM SYSTEM
    η_thermal = wnet / q_in                    #THERMAL EFFECIENCY
    Primary_inlet_temp = θ_1[1]                       #INLET HELIUM TEMPERATURE AT REACTOR HEAT EXCHANGERS
    return BraytonOutput(η_thermal, abs(wnet), abs(wc), abs(wt), abs(wc / wt), abs(q_in), abs(q_out), Primary_inlet_temp)
end

struct circuitOutput{T<:Real}
    T_in::T
    T_out::T
    P_in::T
    P_out::T
    work::T
    heat::T
end

function compressor_stages(par::ParametersActor)
    nstages = par.Nc
    rp = par.rp
    Tin = par.Tmin
    Pin = par.Pmax / par.rp
    effC = 0.89
    cp = 5.1926e3
    cv = 3.1156e3

    a1c = (effC * (rp^(1 / nstages))^(cv / cp) - (rp^(1 / nstages))^(cv / cp) + rp^(1 / nstages)) / (effC * (rp^(1 / nstages))^(cv / cp))
    theta_L = [Tin; Pin]

    A_c = Diagonal([a1c, rp^(1 / nstages)])
    B_c = Diagonal([1 / a1c, 1])
    edot = cp * (a1c - 1) * Tin

    theta_ci = (A_c * B_c)^(nstages - 1) * A_c * theta_L
    T_out = theta_ci[1]
    P_out = theta_ci[2]
    work_in = nstages * edot
    heat_out = -(nstages - 1) * edot
    return circuitOutput(Tin, T_out, Pin, P_out, work_in, heat_out)
end

function turbine_stages(par::ParametersActor)
    nstages = par.Nt
    rp = par.rp
    Tin = par.Tmax
    Pin = par.Pmax

    effT = 0.93
    cp = 5.1926e3
    cv = 3.1156e3
    a1t = (effT * (rp^(1 / nstages))^(cv / cp) - effT * rp^(1 / nstages) + rp^(1 / nstages)) / rp^(1 / nstages)
    theta_H = [Tin; Pin]

    A_t = Diagonal([a1t, rp^(-1 / nstages)])
    B_t = Diagonal([1 / a1t, 1])
    edot = cp * (a1t - 1) * Tin

    theta_ti = (A_t * B_t)^(nstages - 1) * A_t * theta_H
    T_out = theta_ti[1]
    P_out = theta_ti[2]
    work_out = nstages * edot
    heat_in = -(nstages - 1) * edot
    return circuitOutput(Tin, T_out, Pin, P_out, abs(work_out), abs(heat_in))
end

struct ihts_output{T<:Real}
    Tmin::T
    Tin::T
    Tmax::T
    Tout_cycle::T
    pump_work::T
    HX_q::T
    err::T
end

function ihts_heat_exchanger(pump_η::Real, ϵ_hx::Real, Niter::Int64, Tin_cycle::Real, mflow_cycle::Real, cp_cycle::Real, mratio::Real, P_sys::Real, rp_pump::Real, cp_sys::Real, kcoeff_sys::Real)
    mcp_cyc = mflow_cycle * cp_cycle
    mflow_sys = mflow_cycle * mratio
    mcp_sys = mflow_sys * cp_sys

    Cmin = min(mcp_cyc, mcp_sys)
    Tmin_sys = Tin_cycle + 50
    for i in 1:Niter
        oldTmin = Tmin_sys

        Tin_sys = (Tmin_sys) + (Tmin_sys .* rp_pump^kcoeff_sys - Tmin_sys) / pump_η #wall inlet temp

        Tout_sys = Tin_sys + P_sys ./ mcp_sys    #wall outlet temperature

        dT_max = (Tout_sys - Tin_cycle)      #heat exchanger max temperature rise

        qact = ϵ_hx * Cmin .* dT_max
        Tmin_sys = Tout_sys - qact ./ (mcp_sys)
        Tout_cycle = Tin_cycle + qact ./ mcp_cyc
        error_Tmin = abs(Tmin_sys - oldTmin)

        pump_work = mcp_sys .* (Tin_sys - Tmin_sys) / pump_η
        if abs(error_Tmin) < 0.000001 || i == Niter
            return ihts_output(Tmin_sys, Tin_sys, Tout_sys, Tout_cycle, pump_work, qact, error_Tmin)
        end
    end
end

function regenHX(Tci::Real, Thi::Real, eff)
    dT = (Thi - Tci) * eff
    Tco = Tci + dT
    Tho = Thi - dT
    return Tco, Tho
end

function hxeff(Thi::Real, Ch::Real, Tci::Real, Cc::Real, eff::Real)
    Cmin = min(Ch, Cc)
    dt_mx = Thi - Tci
    Q_act = Cmin .* dt_mx .* eff
    dtH = Q_act ./ Ch
    dtC = Q_act ./ Cc
    Tho = Thi - dtH
    Tco = Tci + dtC
    return Q_act, Tho, Tco
end

function getMCP(dd::IMAS.dd)
    cp_he = 5.1926e3

    bop = dd.balance_of_plant
    Tb_max = @ddtime(bop.heat_transfer.breeder.outlet_temperature)
    Tb_min = @ddtime(bop.heat_transfer.breeder.HX_outlet_temperature)

    Tbreeder_ave = (Tb_max + Tb_min) / 2
    cp_pbli, rho_pbli = pbLi_props(Tbreeder_ave)

    mflow_breeder = @ddtime(bop.heat_transfer.breeder.flow_rate)
    mflow_blanket = @ddtime(bop.heat_transfer.wall.flow_rate)
    mflow_divertor = @ddtime(bop.heat_transfer.divertor.flow_rate)
    mflow_cycle = @ddtime(bop.thermal_cycle.flow_rate)

    mcp_breeder = mflow_breeder * cp_pbli
    mcp_blanket = mflow_blanket * cp_he
    mcp_divertor = mflow_divertor * cp_he
    mcp_cycle = mflow_cycle * cp_he
    return mcp_blanket, mcp_divertor, mcp_breeder, mcp_cycle
end

"""
    pbLi_props(temperature::Real)

temperature input in Kelvin
"""
function pbLi_props(temperature::Real)
    specific_heat = (0.195 - 9.116 * 1e-6 .* temperature) .* 1000   #J/kgK
    density = 10520.35 - 1.19051 .* temperature
    return [specific_heat, density]
end

function waste(dd::IMAS.dd, nm::String, sys_power)
    if nm == "wall" || nm == "blk"
        @ddtime(dd.balance_of_plant.heat_transfer.wall.heat_delivered = 0.0)
        @ddtime(dd.balance_of_plant.heat_transfer.wall.heat_waste = sys_power)
    elseif nm == "divertor" || nm == "div"
        @ddtime(dd.balance_of_plant.heat_transfer.divertor.heat_delivered = 0.0)
        @ddtime(dd.balance_of_plant.heat_transfer.divertor.heat_waste = sys_power)
    elseif nm == "breeder"
        @ddtime(dd.balance_of_plant.heat_transfer.breeder.heat_delivered = 0.0)
        @ddtime(dd.balance_of_plant.heat_transfer.breeder.heat_waste = sys_power)
    end
end

function evalBrayton(bout::BraytonOutput, ihts_par::ParametersActor, dd::IMAS.dd)
    cycle_minTemp = bout.T_HX
    Tmax_blk = ihts_par.blanket_max_temp
    Tmax_div = ihts_par.divertor_max_temp
    Tmax_breeder = ihts_par.breeder_hi_temp
    cp = 5.1926e3
    bop = dd.balance_of_plant

    blanket_power = @ddtime(bop.heat_transfer.wall.heat_load) + @ddtime(bop.heat_transfer.wall.circulator_power)
    breeder_power = @ddtime(bop.heat_transfer.breeder.heat_load) + @ddtime(bop.heat_transfer.breeder.circulator_power)
    divertor_power = @ddtime(bop.heat_transfer.divertor.heat_load) + @ddtime(bop.heat_transfer.divertor.circulator_power)

    mcp_blanket, mcp_divertor, mcp_breeder, mcp_cycle = getMCP(dd)

    if cycle_minTemp < Tmax_blk
        Q_act, Tho, Tco = hxeff(Tmax_blk, mcp_blanket, cycle_minTemp, mcp_cycle, ihts_par.blanket_HX_ϵ)

        @ddtime(dd.balance_of_plant.heat_transfer.wall.heat_delivered = Q_act)
        @ddtime(dd.balance_of_plant.heat_transfer.wall.heat_waste = blanket_power - Q_act)

        cycle_minTemp = Tco
    else
        waste(dd, "wall", blanket_power)
    end

    if cycle_minTemp < Tmax_div
        Q_div, Tho, Tco = hxeff(Tmax_div, mcp_divertor, cycle_minTemp, mcp_cycle, ihts_par.divertor_HX_ϵ)

        @ddtime(dd.balance_of_plant.heat_transfer.divertor.heat_delivered = Q_div)
        @ddtime(dd.balance_of_plant.heat_transfer.divertor.heat_waste = divertor_power - Q_div)

        cycle_minTemp = Tco
    else
        waste(dd, "wall", breeder_power)
    end

    if cycle_minTemp < Tmax_breeder
        Q_breed, Tho, Tco = hxeff(Tmax_breeder, mcp_breeder, cycle_minTemp, mcp_cycle, ihts_par.breeder_HX_ϵ)

        @ddtime(dd.balance_of_plant.heat_transfer.breeder.heat_delivered = Q_breed)
        @ddtime(dd.balance_of_plant.heat_transfer.breeder.heat_waste = breeder_power - Q_breed)
        cycle_minTemp = Tco
    else
        waste(dd, "breeder", breeder_power)
    end

    return Tco
end

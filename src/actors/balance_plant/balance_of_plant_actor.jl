#= =================== =#
#  ActorBalanceOfPlant  #
#= =================== =#
Base.@kwdef mutable struct FUSEparameters__ActorBalanceOfPlant{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(Nothing)
    _name::Symbol = :not_set
    generator_conversion_efficiency::Entry{T} = Entry(T, "-", "Efficiency of the generator"; default=0.95) #  Appl. Therm. Eng. 76 (2015) 123–133, https://doi.org/10.1016/j.applthermaleng.2014.10.093
    do_plot::Entry{Bool} = Entry(Bool, "-", "plot"; default=false)
end

mutable struct ActorBalanceOfPlant <: FacilityAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorBalanceOfPlant
    act::ParametersAllActors
    thermal_cycle_actor::ActorThermalCycle
    IHTS_actor::ActorHeatTransfer
    power_needs_actor::ActorPowerNeeds
end

"""
    ActorBalanceOfPlant(dd::IMAS.dd, act::ParametersAllActors; kw...)

Balance of plant actor that estimates the net electrical power output by comparing the balance of plant electrical needs with the electricity generated from the thermal cycle.
!!! note 
    Stores data in `dd.balance_of_plant`
"""
function ActorBalanceOfPlant(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorBalanceOfPlant(kw...)
    actor = ActorBalanceOfPlant(dd, par, act)
    step(actor)
    finalize(actor)
    return actor
end

function ActorBalanceOfPlant(dd::IMAS.dd, par::FUSEparameters__ActorBalanceOfPlant, act::ParametersAllActors; kw...)
    logging_actor_init(ActorBalanceOfPlant)
    par = par(kw...)

    # set the time
    @ddtime(dd.balance_of_plant.time = dd.global_time)

    breeder_hi_temp, breeder_low_temp, cycle_tmax = initial_temperatures(act.ActorThermalCycle.power_cycle_type)

    IHTS_actor = ActorHeatTransfer(dd, act.ActorHeatTransfer, act; breeder_hi_temp, breeder_low_temp)
    thermal_cycle_actor = ActorThermalCycle(dd, act.ActorThermalCycle, act)
    power_needs_actor = ActorPowerNeeds(dd, act.ActorPowerNeeds, act)
    return ActorBalanceOfPlant(dd, par, act, thermal_cycle_actor, IHTS_actor, power_needs_actor)
end


"""
    initial_temperatures(power_cycle_type::Symbol)

    intializes initial temperatures of the coolant loops based on cycle type
"""
function initial_temperatures(power_cycle_type::Symbol)
    breeder_tmax = 1100.0 + 273.15
    breeder_tmin = 550.0 + 273.15
    if power_cycle_type == :rankine_only
        breeder_tmax = 650.0 + 273.15
        breeder_tmin = 185.0 + 273.15
    end
    cycle_tmax = breeder_tmax - 50.0
    return breeder_tmax, breeder_tmin, cycle_tmax
end

function _step(actor::ActorBalanceOfPlant)
    dd = actor.dd
    par = actor.par
    bop = dd.balance_of_plant

    bop_thermal = bop.thermal_cycle
    @ddtime(bop_thermal.generator_conversion_efficiency = par.generator_conversion_efficiency)

    finalize(step(actor.IHTS_actor))
    finalize(step(actor.thermal_cycle_actor))
    finalize(step(actor.power_needs_actor))

    if par.do_plot
        core = sys_coords(dd)
        pl = plot(core)
        regen_xpt = 13.5
        blk_hx_xpoint = 18.5
        div_hx_xpoint = 23.5
        breeder_hx_xpoint = 28.5
        turb_xpt = 33.0

        if dd.balance_of_plant.power_cycle_type == "complex_brayton"
            breeder_hx_xpoint = 26.5
            blk_hx_xpoint = 13.5
            div_hx_xpoint = 18.5
            regen_xpt = 22.5
            turb_xpt = 30.0
        end

        wall_path = wall_cooling_route(core)
        hx1 = init_hx("wall hx", blk_hx_xpoint, -8.0, 2.0, 4.0)
        plot!(hx1)
        wall_path = attach2hx(wall_path, hx1)
        plot!(wall_path.x, wall_path.y, color=:black, linewidth=5, label=nothing)
        plot!(wall_path.x, wall_path.y, color=:steelblue, linewidth=2.5, label=wall_path.name)

        divertor_path = divertor_flow_path(core)
        hx2 = init_hx("divertor hx", div_hx_xpoint, -8.0, 2.0, 4.0)
        plot!(hx2)
        divertor_path = attach2hx(divertor_path, hx2)
        plot!(divertor_path.x, divertor_path.y, color=:black, linewidth=5, label=nothing)
        plot!(divertor_path.x, divertor_path.y, color=:lightpink, linewidth=2.5, label=divertor_path.name, legend=false)

        # breeder_hx_xpoint = 23.5
        # regen_xpt = 28.5
        # turb_xpt = 33.0
        # if dd.balance_of_plant.power_cycle_type=="complex_brayton"
        #     breeder_hx_xpoint = 26.5
        #     regen_xpt = 22.5
        #     turb_xpt = 30
        # end
        breeder_path = breeder_cooling_route(core)
        hx3 = init_hx("breeder hx", breeder_hx_xpoint, -8.0, 2.0, 4.0)
        plot!(hx3)
        breeder_path = attach2hx(breeder_path, hx3)
        plot!(breeder_path.x, breeder_path.y, color=:black, linewidth=5, label=nothing)
        plot!(breeder_path.x, breeder_path.y, color=:red, linewidth=3, label=breeder_path.name, legend=:topleft)
        ylims!(-17, 13)

        regen = init_hx("regen", regen_xpt, -12.0, 2.0, 4.0)
        plot!(regen)

        T = turb("turbine", coords("turbine", [turb_xpt], [-12.0]), 3.0, 4.0)
        C1 = comp("comp", coords("c", [0.0], [-12.0]), 2.0, 2.5)
        C2 = comp("comp", coords("c", [5.0], [-12.0]), 2.0, 2.5)
        C3 = comp("comp", coords("c", [10.0], [-12.0]), 2.0, 2.5)

        ic1 = intercooler("ic", coords("ic", [2.5], [-10.0]), 1.5, 3.5)
        ic2 = intercooler("ic", coords("ic", [7.5], [-10.0]), 1.5, 3.5)
        v1 = [C1, ic1, C2, ic2, C3, hx1, hx2, regen, hx3, T]
        v2 = [C1, ic1, C2, ic2, C3, regen, hx1, hx2, hx3, T]
        cp = v2
        if dd.balance_of_plant.power_cycle_type == "complex_brayton"
            cp = v1
        end
        p = cyclePath(cp)
        plot!(p.x, p.y, color=:black, label="cycle path")
        plot!(T)
        plot!(C1)
        plot!(C2)
        plot!(C3)
        plot!(ic1)
        plot!(ic2)
        plot!(pl, legend=true, grid=false, axis=nothing, border=:none)
        display(pl)
    end

    return actor
end
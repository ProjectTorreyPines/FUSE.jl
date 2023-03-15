#= =================== =#
#  ActorBalanceOfPlant  #
#= =================== =#
Base.@kwdef mutable struct FUSEparameters__ActorBalanceOfPlant{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(Nothing)
    _name::Symbol = :not_set
    needs_model::Switch{Symbol} = Switch(Symbol, [:gasc, :EU_DEMO], "-", "Power plant electrical needs model"; default=:EU_DEMO)
    generator_conversion_efficiency::Entry{T} = Entry(T, "-", "Efficiency of the steam cycle, thermal to electric"; default=0.9)
    do_plot::Entry{Bool} = Entry(Bool, "-", "plot"; default=false)
end

mutable struct ActorBalanceOfPlant <: FacilityAbstractActor
    dd::IMAS.dd
    par::FUSEparameters__ActorBalanceOfPlant
    act::ParametersAllActors
    thermal_cycle_actor::ActorThermalCycle
    IHTS_actor::ActorHeatTransfer
end

"""
    ActorBalanceOfPlant(dd::IMAS.dd, act::ParametersAllActors; kw...)

Balance of plant actor that estimates the net electrical power output by comparing the balance of plant electrical needs with the electricity generated from the thermal cycle.

* `needs_model = :gasc` simply assumes that the power to balance a plant is 7% of the electricity generated.
* `needs_model = :EU_DEMO` subdivides the power plant electrical needs to [:cryostat, :tritium_handling, :pumping] using  EU-DEMO numbers.

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

    breeder_hi_temp, breeder_low_temp, cycle_tmax = ihts_specs(act.ActorThermalCycle.power_cycle_type)

    IHTS_actor = ActorHeatTransfer(dd, act; breeder_hi_temp, breeder_low_temp)
    thermal_cycle_actor = ActorThermalCycle(dd, act; Tmax=cycle_tmax, rp=3.0)
    return ActorBalanceOfPlant(dd, par, act, thermal_cycle_actor, IHTS_actor)
end

function ihts_specs(power_cycle_type::Symbol)
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
    bop_thermal.generator_conversion_efficiency = par.generator_conversion_efficiency .* ones(length(bop.time))
    bop_thermal.power_electric_generated = bop_thermal.net_work .* par.generator_conversion_efficiency .* ones(length(bop.time))

    @ddtime(bop_thermal.total_useful_heat_power = @ddtime(bop.heat_transfer.wall.heat_delivered) + @ddtime(bop.heat_transfer.divertor.heat_delivered) + @ddtime(bop.heat_transfer.breeder.heat_delivered))

    bop_electric = bop.power_electric_plant_operation

    ## heating and current drive systems
    sys = resize!(bop_electric.system, "name" => "H&CD", "index" => 1)
    sys.power = zeros(length(bop.time))
    for (idx, hcd_system) in enumerate(intersect([:nbi, :ec_launchers, :ic_antennas, :lh_antennas], keys(dd)))
        sub_sys = resize!(sys.subsystem, "name" => string(hcd_system), "index" => idx)
        sub_sys.power = electricity(getproperty(dd, hcd_system), bop.time)
        sys.power .+= sub_sys.power
    end

    ## balance of plant systems
    if par.needs_model == :gasc
        sys = resize!(bop_electric.system, "name" => "BOP_gasc", "index" => 2)
        sys.power = 0.07 .* bop_thermal.power_electric_generated

    elseif par.needs_model == :EU_DEMO
        # More realistic DEMO numbers
        bop_systems = [:cryostat, :tritium_handling, :pumping, :pf_active] # index 2 : 5
        for (idx, system) in enumerate(bop_systems)
            sys = resize!(bop_electric.system, "name" => string(system), "index" => (idx + 1))
            sys.power = electricity(system, bop.time)
        end
    else
        error("ActorBalanceOfPlant: par.needs_model = $(par.needs_model) not recognized")
    end


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
        display(pl)
    end

    return actor
end

function heating_and_current_drive_calc(system_unit, time_array::Vector{<:Real})
    power_electric_total = zeros(length(time_array))
    for item_unit in system_unit
        efficiency = prod([getproperty(item_unit.efficiency, i) for i in keys(item_unit.efficiency)])
        power_electric_total .+= IMAS.get_time_array(item_unit.power_launched, :data, time_array, :constant) ./ efficiency
    end
    return power_electric_total
end

function electricity(nbi::IMAS.nbi, time_array::Vector{<:Real})
    return heating_and_current_drive_calc(nbi.unit, time_array)
end

function electricity(ec_launchers::IMAS.ec_launchers, time_array::Vector{<:Real})
    return heating_and_current_drive_calc(ec_launchers.beam, time_array)
end

function electricity(ic_antennas::IMAS.ic_antennas, time_array::Vector{<:Real})
    return heating_and_current_drive_calc(ic_antennas.antenna, time_array)
end

function electricity(lh_antennas::IMAS.lh_antennas, time_array::Vector{<:Real})
    return heating_and_current_drive_calc(lh_antennas.antenna, time_array)
end

function electricity(symbol::Symbol, time_array::Vector{<:Real})
    return electricity(Val{symbol}, time_array)
end

# Dummy functions values taken from DEMO 2017  https://iopscience.iop.org/article/10.1088/0029-5515/57/1/016011
function electricity(::Type{Val{:cryostat}}, time_array::Vector{<:Real})
    return 30e6 .* ones(length(time_array)) # MWe
end

function electricity(::Type{Val{:tritium_handling}}, time_array::Vector{<:Real})
    return 15e6 .* ones(length(time_array)) # MWe
end

function electricity(::Type{Val{:pumping}}, time_array::Vector{<:Real})
    return 80e6 .* ones(length(time_array)) # MWe    (Note this should not be a constant!)
end

function electricity(::Type{Val{:pf_active}}, time_array::Vector{<:Real})
    return 0e6 .* ones(length(time_array)) # MWe    (Note this should not be a constant!)
end

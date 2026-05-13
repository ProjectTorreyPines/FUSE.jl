"""
    init_balance_of_plant!(dd::IMAS.DD, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.DD=IMAS.dd())

Initialize `dd.balance_of_plant` starting from `ini` and `act` parameters.
Accepts any `IMAS.DD` subtype (e.g. IFEdd's `dd_ife`) that embeds `balance_of_plant`.
"""
function init_balance_of_plant!(dd::IMAS.DD, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.DD=IMAS.dd())
    TimerOutputs.reset_timer!("init_bop")
    TimerOutputs.@timeit timer "init_bop" begin
        init_from = ini.general.init_from

        if init_from == :ods
            if !ismissing(dd1.balance_of_plant.power_plant, :power_cycle_type)
                dd.balance_of_plant = deepcopy(dd1.balance_of_plant)
            else
                init_from = :scalars
            end
        end

        if init_from == :scalars
            dd.balance_of_plant.power_plant.power_cycle_type = string(ini.bop.cycle_type)
        end

        return dd
    end
end
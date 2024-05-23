"""
    init_balance_of_plant!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())

Initialize `dd.balance_of_plant` starting from `ini` and `act` parameters
"""
function init_balance_of_plant!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())
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
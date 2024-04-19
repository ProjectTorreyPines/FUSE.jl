"""
    init_balance_of_plant!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())

Initialize `dd.balance_of_plant thermal cycle type from ini.bop`
"""
function init_balance_of_plant!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())
    TimerOutputs.reset_timer!("init_bop")
    TimerOutputs.@timeit timer "init_bop" begin
        init_from = ini.general.init_from
        dd.balance_of_plant.power_plant.power_cycle_type = String(ini.bop.cycle_type)

        return dd
    end
end
# #= ==================== =#
# #  init stability IDS  #
# #= ==================== =#
# """
#     init_stability(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)

# Initialize `dd.stability` starting from `ini` and `act` parameters
# """
# function init_stability(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors)
#     TimerOutputs.reset_timer!("init_stability")
#     TimerOutputs.@timeit timer "init_stability" begin
#         init_from = ini.general.init_from

#         if init_from == :ods
#             dd1 = IMAS.json2imas(ini.ods.filename)
#             if !ismissing(dd1.stability, :time) && length(keys(dd1.stability.time)) > 0
#                 dd.global_time = max(dd.global_time, maximum(dd1.stability.time))
#                 dd.stability = dd1.stability
#                 #stb = dd.stability.time_slice[]
#             else
#                 init_from = :scalars
#             end
#         end
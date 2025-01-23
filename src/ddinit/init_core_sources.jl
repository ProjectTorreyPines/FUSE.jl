"""
    unique_core_sources_names!(cs::IMAS.core_sources)

add trailing number to sources that have the same name
"""
function unique_core_sources_names!(cs::IMAS.core_sources)
    tmp = Dict{String,Int}()
    for source in cs.source
        if source.identifier.name ∉ keys(tmp)
            tmp[source.identifier.name] = 1
        else
            tmp[source.identifier.name] += 1
        end
    end
    for source in cs.source
        if tmp[source.identifier.name] == 1
            tmp[source.identifier.name] = 0
        end
    end
    for source in reverse(cs.source)
        if source.identifier.name ∈ keys(tmp) && tmp[source.identifier.name] > 0
            tmp[source.identifier.name] -= 1
            source.identifier.name = "$(source.identifier.name)_$(tmp[source.identifier.name]+1)"
        end
    end
end

"""
    init_core_sources!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())

Initialize `dd.nbi`, `dd.ec_launchers`, `dd.ic_antennas`, `dd.lh_antennas` starting from `ini` and `act` parameters
"""
function init_core_sources!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())
    TimerOutputs.reset_timer!("init_core_sources")
    TimerOutputs.@timeit timer "init_core_sources" begin
        init_from = ini.general.init_from

        # core_sources
        if init_from == :ods && IMAS.hasdata(dd1.core_sources, :time) && length(dd1.core_sources.time) > 0
            dd.core_sources = deepcopy(dd1.core_sources)
            unique_core_sources_names!(dd.core_sources)
            if isempty(dd1.ec_launchers.beam) && findfirst(:ec, dd.core_sources.source) !== missing
                @assert act.ActorHCD.ec_model == :none "Init using sources from ODS requires `act.ActorHCD.ec_model = :none` or data in `ods.ec_launchers.beam`"
            end
            if isempty(dd1.ic_antennas.antenna) && findfirst(:ic, dd.core_sources.source) !== missing
                @assert act.ActorHCD.ic_model == :none "Init using sources from ODS requires `act.ActorHCD.ic_model = :none` or data in `ods.ic_launchers.antenna`"
            end
            if isempty(dd1.lh_antennas.antenna) && findfirst(:lh, dd.core_sources.source) !== missing
                @assert act.ActorHCD.lh_model == :none "Init using sources from ODS requires `act.ActorHCD.lh_model = :none` or data in `ods.lh_antennas.antenna`"
            end
            if isempty(dd1.nbi.unit) && findfirst(:nbi, dd.core_sources.source) !== missing
                @assert act.ActorHCD.nb_model == :none "Init using sources from ODS requires `act.ActorHCD.nb_model = :none` or data in `ods.nbi.unit`"
            end
            if !isempty(dd.pellets.time_slice) && !isempty(dd.pellets.time_slice[].pellet) && findfirst(:pellet, dd.core_sources.source) !== missing
                @assert act.ActorHCD.pellet_model == :none "Init using sources from ODS requires `act.ActorHCD.pellet_model = :none` or data in `dd.pellets.time_slice[].pellet`"
            end
        else
            empty!(dd.core_sources) # needed for power_scaling_cost_function
        end

        return dd
    end
end


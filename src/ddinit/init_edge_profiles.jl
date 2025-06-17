"""
    init_edge_profiles!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())

Initialize `dd.edge_profiles` starting from `ini` and `act` parameters
"""
function init_edge_profiles!(dd::IMAS.dd{T}, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd()) where {T<:Real}
    TimerOutputs.reset_timer!("init_edge_profiles")
    TimerOutputs.@timeit timer "init_edge_profiles" begin
        init_from = ini.general.init_from

        if init_from == :ods
            if !isempty(dd1.edge_profiles.profiles_1d)
                dd.edge_profiles = deepcopy(dd1.edge_profiles)
            else
                init_from = :scalars
            end
        end

        if init_from == :scalars
            cp1d = dd.core_profiles.profiles_1d[]
            ep1d = resize!(dd.edge_profiles.profiles_1d)

            # For now, we only Initialize values at the separatrix, where ρpol_norm = 1
            ep1d.grid.rho_pol_norm = [one(T)]

            # density
            ep1d.electrons.density = [cp1d.electrons.density[end]]

            # temperature
            ep1d.electrons.temperature = [cp1d.electrons.temperature[end]]
        end

        return dd
    end
end

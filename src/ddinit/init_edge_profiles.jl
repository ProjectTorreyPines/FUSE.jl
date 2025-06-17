"""
    init_edge_profiles!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())

Initialise `dd.edge_profiles` starting from `ini` and `act` parameters
"""
function init_edge_profiles!(dd::IMAS.dd, ini::ParametersAllInits, act::ParametersAllActors, dd1::IMAS.dd=IMAS.dd())
    TimerOutputs.reset_timer!("init_edge_profiles")
    TimerOutputs.@timeit timer "init_edge_profiles" begin
        init_from = ini.general.init_from

        if init_from == :ods
            if !isempty(dd1.edge_profiles.profiles_1d)
                dd.edge_profiles = deepcopy(dd1.edge_profiles)
                ep1d = dd.edge_profiles.profiles_1d[]

            else
                init_from = :scalars
            end

        end

        if !isempty(dd.equilibrium.time_slice)
            equil = dd.equilibrium.time_slice[]
        else
            equil = ini.equilibrium
        end

        if init_from == :scalars
            init_edge_profiles!(
                dd.edge_profiles,
                equil,
                dd.summary;
            )
        end

        return dd
    end
end

function init_edge_profiles!(
    ep::IMAS.edge_profiles,
    equil::Union{IMAS.equilibrium__time_slice,FUSEparameters__equilibrium},
    summary::IMAS.summary;
    )

    ep1d = resize!(ep.profiles_1d)

    # Initialising edge profiles - WIP

    # grid
    # For now, we only initialise values at the separatrix, where ρpol_norm = 1
    ep1d.grid.rho_pol_norm = [1.0]

    # electrons - values should be set in a more meaningful way. For now, these are all just placeholders.

    # density [m^-3]
    ep1d.electrons.density = [1.0e+20]

    # temperature [eV]
    ep1d.electrons.temperature = [500.0]

    return ep
end

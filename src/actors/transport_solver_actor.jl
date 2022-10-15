import NLsolve

#= ========================== =#
#     transport solver actor   #
#= ========================== =#

mutable struct ActorTransportSolver <: PlasmaAbstractActor
    dd::IMAS.dd
    par::ParametersActor
    actor_ct::ActorCoreTransport
end

function ParametersActor(::Type{Val{:ActorTransportSolver}})
    par = ParametersActor(nothing)
    par.evolve_Ti = Switch([:flux_match, :fixed], "", "How to evolve the average ion temperature "; default=:flux_match)
    par.evolve_Te = Switch([:flux_match, :fixed], "", "How to evolve the electron temperature"; default=:flux_match)
    par.evolve_densities = Entry(Union{Dict,Symbol}, "", "OrderedDict to specify which ion species are evolved, kept constant or used for quasi neutarlity"; default=:fixed)
    par.evolve_rotation = Switch([:flux_match, :fixed], "", "How to evolve the electron temperature"; default=:fixed)
    par.rho_transport = Entry(AbstractVector{<:Real}, "", "Rho transport grid"; default=0.2:0.1:0.8)
    par.max_iterations = Entry(Int, "", "Maximum optimizer iterations"; default=50)
    par.step_size = Entry(Real, "", "Step size for each algorithm iteration (note this has a different meaning for each algorithm)"; default=0.25)
    par.optimizer_algorithm = Switch([:anderson, :jacobian_based], "", "Optimizing algorithm used for the flux matching"; default=:anderson)
    par.do_plot = Entry(Bool, "", "Plots the flux matching"; default=false)
    par.verbose = Entry(Bool, "", "Print trace and optimization result"; default=false)
    return par
end

"""
    ActorTransportSolver(dd::IMAS.dd, act::ParametersAllActors; kw...)

The ActorTransportSolver evalutes the transport fluxes and source fluxes and minimizes the flux_match error
"""
function ActorTransportSolver(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorTransportSolver(kw...)
    actor = ActorTransportSolver(dd, par, act)
    step(actor)
    finalize(actor)
    return actor
end

function ActorTransportSolver(dd::IMAS.dd, par::ParametersActor, act::ParametersAllActors; kw...)
    logging_actor_init(ActorTransportSolver)
    par = par(kw...)
    actor_ct = ActorCoreTransport(dd, act.ActorCoreTransport, act; par.rho_transport)
    ActorTransportSolver(dd, par, actor_ct)
end

"""
    _step(actor::ActorTransportSolver)

ActorTransportSolver step
"""
function _step(actor::ActorTransportSolver)
    dd = actor.dd
    par = actor.par

    if par.do_plot
        p = plot(dd.core_profiles, label="before")
    end

    z_init = pack_z_profiles(dd, par)
    old_log_level = log_topics[:actors]
    res = try
        log_topics[:actors] = Logging.Warn
        if par.optimizer_algorithm == :anderson
            res = NLsolve.nlsolve(z -> flux_match_errors(actor, z), z_init, show_trace=par.verbose, method=:anderson, m=5, beta=-par.step_size, iterations=par.max_iterations, ftol=1E-3, xtol=1E-2)
        elseif par.optimizer_algorithm == :jacobian_based
            res = NLsolve.nlsolve(z -> flux_match_errors(actor, z), z_init, show_trace=par.verbose, factor=par.step_size, iterations=par.max_iterations, ftol=1E-3, xtol=1E-2)
        end
        res
    finally
        log_topics[:actors] = old_log_level
    end

    if par.verbose
        display(res)
    end

    flux_match_errors(actor::ActorTransportSolver, res.zero) # res.zero == z_profiles for the smallest error iteration
    if par.evolve_densities !== :fixed
        IMAS.enforce_quasi_neutrality!(dd, [i for (i, evolve) in par.evolve_densities if evolve == :quasi_neutrality][1])
    end

    """
    if !IMAS.is_quasi_neutral(dd)
        IMAS.enforce_quasi_neutrality!(dd,  [i for (i, evolve) in par.evolve_densities if evolve == :quasi_neutrality][1])
    end
    """

    if par.do_plot
        display(plot(dd.core_transport))
        display(plot!(p, dd.core_profiles, label="after"))
    end

    return actor
end

"""
    flux_match_errors(actor::ActorTransportSolver, z_profiles::AbstractVector{<:Float64})

Update the profiles, evaluates neoclassical and turbulent fluxes, sources (ie target fluxes), and returns error between the two
"""
function flux_match_errors(actor::ActorTransportSolver, z_profiles::AbstractVector{<:Float64})
    dd = actor.dd
    par = actor.par

    # modify dd with new z_profiles
    unpack_z_profiles(dd.core_profiles.profiles_1d[], par, z_profiles)

    # evaludate neoclassical + turbulent fluxes
    finalize(step(actor.actor_ct))

    # evaluate sources (ie. target fluxes)
    IMAS.sources!(dd)

    # compare fluxes
    return flux_match_errors(dd, par)
end

function error_transformation!(error, norm)
    error .= error ./ norm
    error .= asinh.(error) # avoid exploding errors, improve robustness
    return error
end

"""
    flux_match_errors(dd::IMAS.dd, par::ParametersActor)

Evaluates the flux_matching errors for the :flux_match species and channels
NOTE: flux matching is done in physical units
"""
function flux_match_errors(dd::IMAS.dd, par::ParametersActor)
    if par.verbose
        flush(stdout)
    end

    cp1d = dd.core_profiles.profiles_1d[]
    total_sources = IMAS.total_sources(dd.core_sources, cp1d)
    total_fluxes = IMAS.total_fluxes(dd.core_transport)

    cs_gridpoints = [argmin((rho_x .- total_sources.grid.rho_tor_norm) .^ 2) for rho_x in par.rho_transport]
    ct_gridpoints = [argmin((rho_x .- total_fluxes.grid_flux.rho_tor_norm) .^ 2) for rho_x in par.rho_transport]
    surface = total_sources.grid.surface[cs_gridpoints]

    error = Float64[]

    if par.evolve_Ti == :flux_match
        norm = 1E4 #[W / m^2]
        err = total_sources.total_ion_power_inside[cs_gridpoints] ./ surface .- total_fluxes.total_ion_energy.flux[ct_gridpoints]
        append!(error, error_transformation!(err, norm))
    end

    if par.evolve_Te == :flux_match
        norm = 1E4 #[W / m^2]
        err = total_sources.electrons.power_inside[cs_gridpoints] ./ surface .- total_fluxes.electrons.energy.flux[ct_gridpoints]
        append!(error, error_transformation!(err, norm))
    end

    if par.evolve_rotation == :flux_match
        norm = 0.1 #[kg / m s^2]
        err = total_sources.torque_tor_inside[cs_gridpoints] ./ surface .- total_fluxes.momentum_tor.flux[ct_gridpoints]
        append!(error, error_transformation!(err, norm))
    end

    if par.evolve_densities != :fixed
        norm = 1E19 #[m^-2 s^-1]
        if par.evolve_densities[:electrons] == :flux_match
            err = total_sources.electrons.particles_inside[cs_gridpoints] ./ surface .- total_fluxes.electrons.particles.flux[ct_gridpoints]
            append!(error, error_transformation!(err, norm))
        end
        for ion in cp1d.ion
            if par.evolve_densities[Symbol(ion.label)] == :flux_match
                error("This is currently not working will fix later")
            end
        end
    end

    return error
end

"""
    pack_z_profiles(dd::IMAS.dd, par::ParametersActor)

Packs the z_profiles based on evolution parameters

NOTE: the order for packing and unpacking is always: [Ti, Te, Rotation, ne, nis...]
"""
function pack_z_profiles(dd::IMAS.dd, par::ParametersActor)
    cp1d = dd.core_profiles.profiles_1d[]
    z_profiles = Float64[]
    cp_gridpoints = [argmin((rho_x .- cp1d.grid.rho_tor_norm) .^ 2) for rho_x in par.rho_transport]

    if par.evolve_Ti == :flux_match
        z_Ti = IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.ion[1].temperature)[cp_gridpoints]
        append!(z_profiles, z_Ti)
    end

    if par.evolve_Te == :flux_match
        z_Te = IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.electrons.temperature)[cp_gridpoints]
        append!(z_profiles, z_Te)
    end

    if par.evolve_rotation == :flux_match
        z_rot = IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.rotation_frequency_tor_sonic)[cp_gridpoints]
        append!(z_profiles, z_rot)
    end

    if par.evolve_densities != :fixed
        check_evolve_densities(dd, par.evolve_densities)
        if par.evolve_densities[:electrons] == :flux_match
            z_ne = IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal)[cp_gridpoints]
            append!(z_profiles, z_ne)
        end
        for ion in cp1d.ion
            if par.evolve_densities[Symbol(ion.label)] == :flux_match
                z_ni = IMAS.calcz(cp1d.grid.rho_tor_norm, ion.density_thermal)[cp_gridpoints]
                append!(z_profiles, z_ni)
            end
        end
    end

    return z_profiles
end

"""
    unpack_z_profiles(cp1d::IMAS.core_profiles__profiles_1d, par::ParametersActor, z_profiles::AbstractVector{<:Real})

Unpacks z_profiles based on evolution parameters

NOTE: The order for packing and unpacking is always: [Ti, Te, Rotation, ne, nis...]
"""
function unpack_z_profiles(cp1d::IMAS.core_profiles__profiles_1d, par::ParametersActor, z_profiles::AbstractVector{<:Real})
    rho_transport = par.rho_transport
    counter = 0
    N = length(rho_transport)

    if par.evolve_Ti == :flux_match
        Ti_new = IMAS.profile_from_z_transport(cp1d.ion[1].temperature, cp1d.grid.rho_tor_norm, rho_transport, z_profiles[counter+1:counter+N])
        counter += N
        for ion in cp1d.ion
            ion.temperature = Ti_new
        end
    end

    if par.evolve_Te == :flux_match
        cp1d.electrons.temperature = IMAS.profile_from_z_transport(cp1d.electrons.temperature, cp1d.grid.rho_tor_norm, rho_transport, z_profiles[counter+1:counter+N])
        counter += N
    end

    if par.evolve_rotation == :flux_match
        cp1d.rotation_frequency_tor_sonic = IMAS.profile_from_z_transport(cp1d.rotation_frequency_tor_sonic .+ 1, cp1d.grid.rho_tor_norm, rho_transport, z_profiles[counter+1:counter+N])
        counter += N
    end

    if par.evolve_densities != :fixed
        if par.evolve_densities[:electrons] == :flux_match
            cp1d.electrons.density_thermal = IMAS.profile_from_z_transport(cp1d.electrons.density_thermal, cp1d.grid.rho_tor_norm, rho_transport, z_profiles[counter+1:counter+N])
            counter += N
        end
        z_ne = IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal)[1:argmin(abs.(cp1d.grid.rho_tor_norm .- rho_transport[end]))]
        for ion in cp1d.ion
            if par.evolve_densities[Symbol(ion.label)] == :flux_match
                ion.density_thermal = IMAS.profile_from_z_transport(ion.density_thermal, cp1d.grid.rho_tor_norm, rho_transport, z_profiles[counter+1:counter+N])
                counter += N
            elseif par.evolve_densities[Symbol(ion.label)] == :match_ne_scale
                ion.density_thermal = IMAS.profile_from_z_transport(ion.density_thermal, cp1d.grid.rho_tor_norm, rho_transport, z_ne)

            end
        end
    end

    # Ensure Quasi neutrality
    if !IMAS.is_quasi_neutral(cp1d) && par.evolve_densities != :fixed
        q_specie = [i for (i, evolve) in par.evolve_densities if evolve == :quasi_neutrality]
        @assert q_specie != Symbol[] "no quasi neutrality specie while quasi neutrality is broken"
        IMAS.enforce_quasi_neutrality!(cp1d, q_specie[1])
    end

    return cp1d
end

"""
    check_evolve_densities(dd::IMAS.dd, evolve_densities::Dict)
"""
function check_evolve_densities(dd::IMAS.dd, evolve_densities::Dict)
    return check_evolve_densities(dd.core_profiles.profiles_1d[], evolve_densities)
end

"""
    check_evolve_densities(cp1d::IMAS.core_profiles__profiles_1d, evolve_densities::Dict)

Checks if the evolve_densities dictionary makes sense and return sensible errors if this is not the case
"""
function check_evolve_densities(cp1d::IMAS.core_profiles__profiles_1d, evolve_densities::Dict)
    dd_species = vcat([Symbol(ion.label) for ion in cp1d.ion], :electrons)
    dd_species = vcat(dd_species, [Symbol(String(ion.label) * "_fast") for ion in cp1d.ion if sum(ion.density_fast) > 0.0])
    # Check if evolve_densities contains all of dd species
    @assert sort([i for (i, evolve) in evolve_densities]) == sort(dd_species) "Not all species are accounted for in the evolve_densities dict : $(sort([i for (i,j) in evolve_densities])) , dd_species : $(sort(dd_species)) ,"
    # Check if there is 1 quasi_neutality specie
    @assert length([i for (i, evolve) in evolve_densities if evolve == :quasi_neutality]) < 2 "Only one specie can be used for quasi neutality matching (or 0 when everything matches ne_scale)"
    # Check if there is a source for the flux_match specie(s)
    # isnt' there yet 
end

"""
    setup_density_evolution_electron_flux_match_rest_ne_scale(dd::IMAS.dd)

Sets up the evolve_density dict to evolve only ne and keep the rest matching the ne_scale lengths
"""
function setup_density_evolution_electron_flux_match_rest_ne_scale(dd::IMAS.dd)
    dd_thermal = [Symbol(ion.label) for ion in dd.core_profiles.profiles_1d[].ion if sum(ion.density_thermal) > 0.0]
    dd_fast = [Symbol(String(ion.label) * "_fast") for ion in dd.core_profiles.profiles_1d[].ion if sum(ion.density_fast) > 0.0]
    return evolve_densities_dict_creation([:electrons], dd_fast, dd_thermal)
end

"""
    evolve_densities_dict_creation(flux_match_species::Vector, fixed_species::Vector, match_ne_scale_species::Vector; quasi_neutality_specie::Union{Symbol,Bool}=false)

Create the density_evolution dict based on input vectors: flux_match_species, fixed_species, match_ne_scale_species, quasi_neutality_specie
"""
function evolve_densities_dict_creation(flux_match_species::Vector, fixed_species::Vector, match_ne_scale_species::Vector; quasi_neutality_specie::Union{Symbol,Bool}=false)
    parse_list = vcat([[sym, :flux_match] for sym in flux_match_species], [[sym, :match_ne_scale] for sym in match_ne_scale_species], [[sym, :fixed] for sym in fixed_species])
    if isa(quasi_neutality_specie, Symbol)
        parse_list = vcat(parse_list, [[quasi_neutality_specie, :quasi_neutality]])
    end
    return Dict(sym => evolve for (sym, evolve) in parse_list)
end

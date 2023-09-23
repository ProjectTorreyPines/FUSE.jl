import NLsolve
using LinearAlgebra
#= ========================== =#
#     transport solver actor   #
#= ========================== =#
Base.@kwdef mutable struct FUSEparameters__ActorFluxMatcher{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    evolve_Ti::Switch{Symbol} = Switch{Symbol}([:flux_match, :fixed], "-", "Evolve ion temperature "; default=:flux_match)
    evolve_Te::Switch{Symbol} = Switch{Symbol}([:flux_match, :fixed], "-", "Evolve electron temperature"; default=:flux_match)
    evolve_densities::Entry{Union{AbstractDict,Symbol}} =
        Entry{Union{AbstractDict,Symbol}}("-", "Dict to specify which ion species are evolved, kept constant, or used to enforce quasi neutarlity"; default=:fixed)
    evolve_rotation::Switch{Symbol} = Switch{Symbol}([:flux_match, :fixed], "-", "Evolve the electron temperature"; default=:fixed)
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "Rho transport grid"; default=0.2:0.1:0.8)
    evolve_pedestal::Entry{Bool} = Entry{Bool}("-", "Evolve the pedestal inside the transport solver"; default=true)
    max_iterations::Entry{Int} = Entry{Int}("-", "Maximum optimizer iterations"; default=200)
    optimizer_algorithm::Switch{Symbol} = Switch{Symbol}([:anderson, :jacobian_based], "-", "Optimizing algorithm used for the flux matching"; default=:anderson)
    step_size::Entry{T} = Entry{T}("-", "Step size for each algorithm iteration (note this has a different meaning for each algorithm)"; default=1.0)
    do_plot::Entry{Bool} = Entry{Bool}("-", "Plots the flux matching"; default=false)
    verbose::Entry{Bool} = Entry{Bool}("-", "Print trace and optimization result"; default=false)
end

mutable struct ActorFluxMatcher{D,P} <: PlasmaAbstractActor
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorFluxMatcher{P}
    actor_ct::ActorFluxCalculator{D,P}
    actor_ped::ActorPedestal{D,P}
end

"""
    ActorFluxMatcher(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evalutes the transport fluxes and source fluxes and minimizes the flux_match error
"""
function ActorFluxMatcher(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorFluxMatcher(dd, act.ActorFluxMatcher, act; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorFluxMatcher(dd::IMAS.dd, par::FUSEparameters__ActorFluxMatcher, act::ParametersAllActors; kw...)
    logging_actor_init(ActorFluxMatcher)
    par = par(kw...)
    actor_ct = ActorFluxCalculator(dd, act.ActorFluxCalculator, act; par.rho_transport)
    actor_ped = ActorPedestal(dd, act.ActorPedestal; ip_from=:equilibrium, βn_from=:core_profiles)
    return ActorFluxMatcher(dd, par, actor_ct, actor_ped)
end

"""
    _step(actor::ActorFluxMatcher)

ActorFluxMatcher step
"""
function _step(actor::ActorFluxMatcher)
    dd = actor.dd
    par = actor.par

    if par.do_plot
        cp1d_before = deepcopy(dd.core_profiles.profiles_1d[])
    end

    z_init = pack_z_profiles(dd, par) .* 100
    old_log_level = log_topics[:actors]
    err_history = Float64[]
    z_history = Vector{Float64}[]
    res = try
        log_topics[:actors] = Logging.Warn
        if par.optimizer_algorithm == :anderson
            res = NLsolve.nlsolve(
                z -> flux_match_errors(actor, z; z_history, err_history),
                z_init;
                show_trace=par.verbose,
                store_trace=par.verbose,
                method=:anderson,
                m=0,
                beta=-par.step_size,
                iterations=par.max_iterations,
                ftol=1E-3,
                xtol=1E-2
            )
        elseif par.optimizer_algorithm == :jacobian_based
            res = NLsolve.nlsolve(
                z -> flux_match_errors(actor, z; z_history, err_history),
                z_init;
                factor=1e-2,
                show_trace=par.verbose,
                store_trace=par.verbose,
                iterations=par.max_iterations,
                ftol=1E-3
            )
        end
        res
    finally
        log_topics[:actors] = old_log_level
    end

    if par.verbose
        display(res)
        parse_and_plot_error(string(res.trace.states))
    end

    flux_match_errors(actor::ActorFluxMatcher, z_history[argmin(err_history)]) # res.zero == z_profiles for the smallest error iteration
    if par.evolve_densities !== :fixed
        IMAS.enforce_quasi_neutrality!(dd, [i for (i, evolve) in par.evolve_densities if evolve == :quasi_neutrality][1])
    end

    """
    if !IMAS.is_quasi_neutral(dd)
        IMAS.enforce_quasi_neutrality!(dd,  [i for (i, evolve) in par.evolve_densities if evolve == :quasi_neutrality][1])
    end
    """

    if par.do_plot
        cp1d = dd.core_profiles.profiles_1d[]
        N_channels = Int(length(z_init) / length(par.rho_transport))
        p = plot(; layout=(N_channels, 2), size=(1000, 1000), xguidefontsize=15, yguidefontsize=14, legendfontsize=12,
            tickfont=font(12, "Computer Modern"), fontfamily="Computer Modern")

        titels = ["Electron temperature", "Ion temperature", "Electron density", "Rotation frequency tor sonic"]
        to_plot_after = [(cp1d.electrons, :temperature), (cp1d.ion[1], :temperature), (cp1d.electrons, :density_thermal), (cp1d, :rotation_frequency_tor_sonic)]
        to_plot_before =
            [(cp1d_before.electrons, :temperature), (cp1d_before.ion[1], :temperature), (cp1d_before.electrons, :density_thermal), (cp1d_before, :rotation_frequency_tor_sonic)]

        for sub in 1:N_channels
            plot!(dd.core_transport; only=sub, subplot=2 * sub - 1, aspect=:equal)
            plot!(to_plot_before[sub][1], to_plot_before[sub][2]; subplot=2 * sub, label="before", linestyle=:dash, color=:black)
            plot!(to_plot_after[sub][1], to_plot_after[sub][2]; subplot=2 * sub, label="after", title=titels[sub], aspect=:equal)
        end
        display(p)
    end

    return actor
end

"""
    flux_match_errors(actor::ActorFluxMatcher, z_profiles::AbstractVector{<:Real})

Update the profiles, evaluates neoclassical and turbulent fluxes, sources (ie target fluxes), and returns error between the two
"""
function flux_match_errors(
    actor::ActorFluxMatcher,
    z_profiles::AbstractVector{<:Real};
    z_history::Vector{Vector{Float64}}=Vector{Float64}[],
    err_history::Vector{Float64}=Float64[]
)
    push!(z_history, z_profiles)
    z_profiles = z_profiles ./ 100
    dd = actor.dd
    par = actor.par

    # evolve pedestal
    if par.evolve_pedestal
        # modify dd with new z_profiles
        actor.actor_ped.par.βn_from = :equilibrium
        _finalize(_step(actor.actor_ped))
        unpack_z_profiles(dd.core_profiles.profiles_1d[], par, z_profiles)
        actor.actor_ped.par.βn_from = :core_profiles
        _finalize(_step(actor.actor_ped))
        unpack_z_profiles(dd.core_profiles.profiles_1d[], par, z_profiles)
    else
        # modify dd with new z_profiles
        unpack_z_profiles(dd.core_profiles.profiles_1d[], par, z_profiles)
    end
    # evaluate sources (ie. target fluxes)
    IMAS.sources!(dd)

    # evaludate neoclassical + turbulent fluxes
    _finalize(_step(actor.actor_ct))

    # compare fluxes
    errors = flux_match_errors(dd, par)
    push!(err_history, norm(errors))
    return errors
end

function error_transformation!(target::T, output::T, norm::Real) where {T<:AbstractVector{<:Real}}
    error = (target .- output) ./ norm
    return asinh.(error)
end

"""
    flux_match_errors(dd::IMAS.dd, par::FUSEparameters__ActorFluxMatcher)

Evaluates the flux_matching errors for the :flux_match species and channels
NOTE: flux matching is done in physical units
"""
function flux_match_errors(dd::IMAS.dd, par::FUSEparameters__ActorFluxMatcher)
    cp1d = dd.core_profiles.profiles_1d[]
    total_sources = IMAS.total_sources(dd.core_sources, cp1d; fields=[:total_ion_power_inside, :power_inside, :particles_inside, :torque_tor_inside])
    total_fluxes = IMAS.total_fluxes(dd.core_transport)

    cs_gridpoints = [argmin((rho_x .- total_sources.grid.rho_tor_norm) .^ 2) for rho_x in par.rho_transport]
    cf_gridpoints = [argmin((rho_x .- total_fluxes.grid_flux.rho_tor_norm) .^ 2) for rho_x in par.rho_transport]

    error = Float64[]

    if par.evolve_Ti == :flux_match
        norm = 1E4 #[W / m^2]
        target = total_sources.total_ion_power_inside[cs_gridpoints] ./ total_sources.grid.surface[cs_gridpoints]
        output = total_fluxes.total_ion_energy.flux[cf_gridpoints]
        append!(error, error_transformation!(target, output, norm))
    end

    if par.evolve_Te == :flux_match
        norm = 1E4 #[W / m^2]
        target = total_sources.electrons.power_inside[cs_gridpoints] ./ total_sources.grid.surface[cs_gridpoints]
        output = total_fluxes.electrons.energy.flux[cf_gridpoints]
        append!(error, error_transformation!(target, output, norm))
    end

    if par.evolve_rotation == :flux_match
        norm = 1E-3 #[kg / m s^2]
        target = total_sources.torque_tor_inside[cs_gridpoints] ./ total_sources.grid.surface[cs_gridpoints]
        output = total_fluxes.momentum_tor.flux[cf_gridpoints]
        append!(error, error_transformation!(target, output, norm))
    end

    if par.evolve_densities != :fixed
        norm = 1E19 #[m^-2 s^-1]
        if par.evolve_densities[:electrons] == :flux_match
            target = total_sources.electrons.particles_inside[cs_gridpoints] ./ total_sources.grid.surface[cs_gridpoints]
            output = total_fluxes.electrons.particles.flux[cf_gridpoints]
            append!(error, error_transformation!(target, output, norm))
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
    pack_z_profiles(dd::IMAS.dd, par::FUSEparameters__ActorFluxMatcher)

Packs the z_profiles based on evolution parameters

NOTE: the order for packing and unpacking is always: [Ti, Te, Rotation, ne, nis...]
"""
function pack_z_profiles(dd::IMAS.dd, par::FUSEparameters__ActorFluxMatcher)
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
    unpack_z_profiles(cp1d::IMAS.core_profiles__profiles_1d, par::FUSEparameters__ActorFluxMatcher, z_profiles::AbstractVector{<:Real})

Unpacks z_profiles based on evolution parameters

NOTE: The order for packing and unpacking is always: [Ti, Te, Rotation, ne, nis...]
"""
function unpack_z_profiles(cp1d::IMAS.core_profiles__profiles_1d, par::FUSEparameters__ActorFluxMatcher, z_profiles::AbstractVector{<:Real})
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
        cp1d.rotation_frequency_tor_sonic =
            IMAS.profile_from_z_transport(cp1d.rotation_frequency_tor_sonic .+ 1, cp1d.grid.rho_tor_norm, rho_transport, z_profiles[counter+1:counter+N])
        counter += N
    end

    if par.evolve_densities != :fixed
        if par.evolve_densities[:electrons] == :flux_match
            cp1d.electrons.density_thermal = IMAS.profile_from_z_transport(cp1d.electrons.density_thermal, cp1d.grid.rho_tor_norm, rho_transport, z_profiles[counter+1:counter+N])
            counter += N
        end
        z_ne = IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal)
        for ion in cp1d.ion
            if par.evolve_densities[Symbol(ion.label)] == :flux_match
                ion.density_thermal = IMAS.profile_from_z_transport(ion.density_thermal, cp1d.grid.rho_tor_norm, rho_transport, z_profiles[counter+1:counter+N])
                counter += N
            elseif par.evolve_densities[Symbol(ion.label)] == :match_ne_scale
                ion.density_thermal = IMAS.profile_from_z_transport(ion.density_thermal, cp1d.grid.rho_tor_norm, cp1d.grid.rho_tor_norm, z_ne)
            end
        end
    end
    # Ensure Quasi neutrality if you are evolving densities
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
function check_evolve_densities(dd::IMAS.dd, evolve_densities::AbstractDict)
    return check_evolve_densities(dd.core_profiles.profiles_1d[], evolve_densities)
end

"""
    check_evolve_densities(cp1d::IMAS.core_profiles__profiles_1d, evolve_densities::Dict)

Checks if the evolve_densities dictionary makes sense and return sensible errors if this is not the case
"""
function check_evolve_densities(cp1d::IMAS.core_profiles__profiles_1d, evolve_densities::AbstractDict)
    dd_species = vcat([Symbol(ion.label) for ion in cp1d.ion], :electrons)
    dd_species = vcat(dd_species, [Symbol(String(ion.label) * "_fast") for ion in cp1d.ion if sum(ion.density_fast) > 0.0])
    # Check if evolve_densities contains all of dd species
    @assert sort([i for (i, evolve) in evolve_densities]) == sort(dd_species) "Not all species are accounted for in the evolve_densities dict : $(sort([i for (i,j) in evolve_densities])) , dd_species : $(sort(dd_species)) ,"
    # Check if there is 1 quasi_neutrality specie
    @assert length([i for (i, evolve) in evolve_densities if evolve == :quasi_neutrality]) < 2 "Only one specie can be used for quasi neutality matching (or 0 when everything matches ne_scale)"
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
    quasi_neutrality_specie = :D
    if :DT ∈ dd_thermal
        quasi_neutrality_specie = :DT
    end
    return evolve_densities_dict_creation([:electrons], dd_fast, dd_thermal; quasi_neutrality_specie)
end

"""
    evolve_densities_dict_creation(flux_match_species::Vector, fixed_species::Vector, match_ne_scale_species::Vector; quasi_neutrality_specie::Union{Symbol,Bool}=false)

Create the density_evolution dict based on input vectors: flux_match_species, fixed_species, match_ne_scale_species, quasi_neutrality_specie
"""
function evolve_densities_dict_creation(flux_match_species::Vector, fixed_species::Vector, match_ne_scale_species::Vector; quasi_neutrality_specie::Union{Symbol,Bool}=false)
    parse_list = vcat([[sym, :flux_match] for sym in flux_match_species], [[sym, :match_ne_scale] for sym in match_ne_scale_species], [[sym, :fixed] for sym in fixed_species])
    if isa(quasi_neutrality_specie, Symbol)
        parse_list = vcat(parse_list, [[quasi_neutrality_specie, :quasi_neutrality]])
    end
    return Dict(sym => evolve for (sym, evolve) in parse_list)
end

function parse_and_plot_error(data::String)
    data = split(data, "\n")[2:end-1]
    array = zeros(length(data))
    for (idx, line) in enumerate(data)
        filtered_arr = filter(x -> !occursin(r"^\s*$", x), split(line, " "))
        array[idx] = parse(Float64, filtered_arr[3])
    end
    return display(plot(array; yscale=:log10, ylabel="log of convergence errror", xlabel="iterations", label=@sprintf("Minimum error =  %.3e ", (minimum(array))), ylim=[1e-4, 10]))
end
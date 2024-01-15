import NLsolve
using LinearAlgebra
import TJLF: InputTJLF

#= ================ =#
#  ActorFluxMatcher  #
#= ================ =#
Base.@kwdef mutable struct FUSEparameters__ActorFluxMatcher{T} <: ParametersActor where {T<:Real}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "ρ transport grid"; default=0.25:0.1:0.85)
    evolve_Ti::Switch{Symbol} = Switch{Symbol}([:flux_match, :fixed], "-", "Ion temperature `:flux_match` or keep `:fixed`"; default=:flux_match)
    evolve_Te::Switch{Symbol} = Switch{Symbol}([:flux_match, :fixed], "-", "Electron temperature `:flux_match` or keep `:fixed`"; default=:flux_match)
    evolve_densities::Entry{Union{AbstractDict,Symbol}} =
        Entry{Union{AbstractDict,Symbol}}(
            "-",
            "Densities `:fixed`, or electron flux-match and rest match ne scale `:flux_match`, or Dict to specify which species are `:flux_match`, kept `:fixed`, used to enforce `:quasi_neutrality`, or scaled to `:match_ne_scale`";
            default=:fixed
        )
    evolve_rotation::Switch{Symbol} = Switch{Symbol}([:flux_match, :fixed], "-", "Rotation `:flux_match` or keep `:fixed`"; default=:fixed)
    evolve_pedestal::Entry{Bool} = Entry{Bool}("-", "Evolve the pedestal inside the transport solver"; default=false)
    optimize_q::Entry{Bool} = Entry{Bool}("-", "Optimize q profile to achieve goal"; default=false)
    find_widths::Entry{Bool} = Entry{Bool}("-", "Runs Turbulent transport actor TJLF finding widths after first iteration"; default=true)
    max_iterations::Entry{Int} = Entry{Int}("-", "Maximum optimizer iterations"; default=300)
    optimizer_algorithm::Switch{Symbol} = Switch{Symbol}([:anderson, :newton, :trust_region], "-", "Optimizing algorithm used for the flux matching"; default=:anderson)
    step_size::Entry{T} = Entry{T}("-", "Step size for each algorithm iteration (note this has a different meaning for each algorithm)"; default=1.0)
    do_plot::Entry{Bool} = Entry{Bool}("-", "Plots the flux matching"; default=false)
    verbose::Entry{Bool} = Entry{Bool}("-", "Print trace and optimization result"; default=false)
end

mutable struct ActorFluxMatcher{D,P} <: PlasmaAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::FUSEparameters__ActorFluxMatcher{P}
    actor_ct::ActorFluxCalculator{D,P}
    actor_ped::ActorPedestal{D,P}
    norms::Vector{Float64}
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
    actor_ped = ActorPedestal(dd, act.ActorPedestal; ip_from=:equilibrium, βn_from=:core_profiles, rho_nml=par.rho_transport[end-1], rho_ped=par.rho_transport[end])
    return ActorFluxMatcher(dd, par, actor_ct, actor_ped, Float64[])
end

"""
    _step(actor::ActorFluxMatcher)

ActorFluxMatcher step
"""
function _step(actor::ActorFluxMatcher)
    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]

    @assert nand(typeof(actor.actor_ct.actor_neoc) <: ActorNoOperation, typeof(actor.actor_ct.actor_turb) <: ActorNoOperation) "Unable to fluxmatch when all transport actors are turned off"

    if par.do_plot
        cp1d_before = deepcopy(dd.core_profiles.profiles_1d[])
    end

    # normalizations for each of the channels
    # updated with exponential forget for subsequent step() calls
    if isempty(actor.norms)
        actor.norms = flux_match_norms(dd, par)
    else
        actor.norms = actor.norms .* 0.5 .+ flux_match_norms(dd, par) .* 0.5
    end

    z_init = scale_z_profiles(pack_z_profiles(cp1d, par)) # scale z_profiles to get smaller stepping

    if par.optimizer_algorithm == :newton
        opts = Dict(:method => :newton, :factor => par.step_size)
    elseif par.optimizer_algorithm == :anderson
        @assert 0.0 < par.step_size <= 1.0
        opts = Dict(:method => :anderson, :m => 5, :beta => -par.step_size)
    elseif par.optimizer_algorithm == :trust_region
        opts = Dict(:method => :trust_region, :factor => par.step_size, :autoscale => true)
    end

    prog = ProgressMeter.ProgressUnknown(; desc="Calls:", enabled=par.verbose)

    z_scaled_history = Vector{NTuple{length(z_init),Float64}}()
    err_history = Float64[]
    old_logging = actor_logging(dd, false)
    res = try
        res = NLsolve.nlsolve(
            z -> flux_match_errors(actor, z; z_scaled_history, err_history, prog),
            z_init;
            show_trace=false,
            store_trace=true,
            extended_trace=false,
            iterations=par.max_iterations,
            ftol=1E-3,
            xtol=1E-2,
            opts...
        )
        flux_match_errors(actor, res.zero) # z_profiles for the smallest error iteration
    finally
        # @show flux_match_errors(actor, z_scaled_history[end])
        actor_logging(dd, old_logging)
    end

    if par.do_plot
        display(plot(err_history; yscale=:log10, ylabel="Log₁₀ of convergence errror", xlabel="Iterations", label=@sprintf("Minimum error =  %.3e ", (minimum(err_history)))))
        display(plot(transpose(hcat(map(z -> collect(unscale_z_profiles(z)), z_scaled_history)...)); xlabel="Iterations", label=""))

        N_channels = Int(floor(length(z_init) / length(par.rho_transport)))
        p = plot(; layout=(N_channels, 2), size=(1000, 1000))

        titles = ["Electron temperature", "Ion temperature", "Electron density", "Rotation frequency tor sonic"]
        to_plot_after = [(cp1d.electrons, :temperature), (cp1d.ion[1], :temperature), (cp1d.electrons, :density_thermal), (cp1d, :rotation_frequency_tor_sonic)]
        to_plot_before =
            [(cp1d_before.electrons, :temperature), (cp1d_before.ion[1], :temperature), (cp1d_before.electrons, :density_thermal), (cp1d_before, :rotation_frequency_tor_sonic)]

        for sub in 1:N_channels
            plot!(dd.core_transport; only=sub, subplot=2 * sub - 1, aspect=:equal)
            plot!(to_plot_before[sub][1], to_plot_before[sub][2]; subplot=2 * sub, label="before", linestyle=:dash, color=:black)
            plot!(to_plot_after[sub][1], to_plot_after[sub][2]; subplot=2 * sub, label="after", title=titles[sub], aspect=:equal)
        end

        display(p)
    end

    # for completely collapsed cases we don't want it to crash in the optimizer so
    cp1d = dd.core_profiles.profiles_1d[]
    if cp1d.electrons.temperature[1] < cp1d.electrons.temperature[end] && par.evolve_Te
        @warn "Profiles completely collpased due to insufficient source versus turbulence"
        te = cp1d.electrons.temperature
        teped = @ddtime(dd.summary.local.pedestal.t_e.value)
        lowest_profile = IMAS.Hmode_profiles(te[end], teped, teped * 1.5, length(te), 1.1, 1.1, 2 * (1 - @ddtime(dd.summary.local.pedestal.position.rho_tor_norm)))
        cp1d.electrons.temperature = lowest_profile
        if par.evolve_Ti
            for ion in cp1d.ion
                ion.temperature = lowest_profile
            end
        end
    end

    return actor
end

"""
    scale_z_profiles(z_profiles)

scale z_profiles to get smaller stepping in nlsolve
"""
function scale_z_profiles(z_profiles)
    return z_profiles .* 100
end

function unscale_z_profiles(z_profiles_scaled)
    return z_profiles_scaled ./ 100
end

function flux_match_errors(
    actor::ActorFluxMatcher,
    z_profiles_scaled::Tuple;
    z_scaled_history::Vector=[],
    err_history::Vector{Float64}=Float64[],
    prog::Any=nothing)
    return flux_match_errors(actor, collect(z_profiles_scaled); z_scaled_history, err_history, prog)
end

"""
    flux_match_errors(actor::ActorFluxMatcher, z_profiles::AbstractVector{<:Real})

Update the profiles, evaluates neoclassical and turbulent fluxes, sources (ie target fluxes), and returns error between the two
"""
function flux_match_errors(
    actor::ActorFluxMatcher,
    z_profiles_scaled::Vector{<:Real};
    z_scaled_history::Vector=[],
    err_history::Vector{Float64}=Float64[],
    prog::Any=nothing)

    dd = actor.dd
    par = actor.par
    cp1d = dd.core_profiles.profiles_1d[]

    # unscale z_profiles
    push!(z_scaled_history, Tuple(z_profiles_scaled))
    z_profiles = unscale_z_profiles(z_profiles_scaled)

    # evolve pedestal
    if par.evolve_pedestal
        unpack_z_profiles(cp1d, par, z_profiles)
        actor.actor_ped.par.βn_from = :core_profiles
        finalize(step(actor.actor_ped))
    end

    # modify dd with new z_profiles
    unpack_z_profiles(cp1d, par, z_profiles)

    # evaluate sources (ie. target fluxes)
    IMAS.sources!(dd)

    if !par.find_widths && length(err_history) > 0 && actor.actor_ct.actor_turb.par.model == :TJLF
        for input_tglf in actor.actor_ct.actor_turb.input_tglfs
            input_tglf.FIND_WIDTH = false
        end
    end
    # evaludate neoclassical + turbulent fluxes
    finalize(step(actor.actor_ct))

    # compare fluxes
    errors = flux_match_errors(dd, par, actor.norms, prog)

    # update error history
    push!(err_history, norm(errors))

    return errors
end

function error_transformation(target::T, output::T, norm::Real) where {T}
    return (target .- output) ./ norm
end

function target_transformation(target::Vector{Float64})
    return sum(abs.(target)) / length(target)
end

"""
    flux_match_norms(dd::IMAS.dd, par::FUSEparameters__ActorFluxMatcher)

Returns vector of normalizations for the different channels that are evolved
"""
function flux_match_norms(dd::IMAS.dd, par::FUSEparameters__ActorFluxMatcher)
    cp1d = dd.core_profiles.profiles_1d[]
    total_sources = IMAS.total_sources(dd.core_sources, cp1d; fields=[:total_ion_power_inside, :power_inside, :particles_inside, :torque_tor_inside])
    cs_gridpoints = [argmin((rho_x .- total_sources.grid.rho_tor_norm) .^ 2) for rho_x in par.rho_transport]

    norms = Float64[]

    if par.evolve_Ti == :flux_match #[W / m^2]
        push!(norms, target_transformation(total_sources.total_ion_power_inside[cs_gridpoints] ./ total_sources.grid.surface[cs_gridpoints]))
    end

    if par.evolve_Te == :flux_match #[W / m^2]
        push!(norms, target_transformation(total_sources.electrons.power_inside[cs_gridpoints] ./ total_sources.grid.surface[cs_gridpoints]))
    end

    if par.evolve_rotation == :flux_match #[kg / m s^2]
        push!(norms, target_transformation(total_sources.torque_tor_inside[cs_gridpoints] ./ total_sources.grid.surface[cs_gridpoints]))
    end

    evolve_densities = evolve_densities_dictionary(cp1d, par)
    if evolve_densities[:electrons] == :flux_match #[m^-2 s^-1]
        push!(norms, target_transformation(total_sources.electrons.particles_inside[cs_gridpoints] ./ total_sources.grid.surface[cs_gridpoints]))
    end

    return norms
end

"""
    flux_match_errors(dd::IMAS.dd, par::FUSEparameters__ActorFluxMatcher, norms::Vector{Float64}, prog::Any)

Evaluates the flux_matching errors for the :flux_match species and channels

NOTE: flux matching is done in physical units
"""
function flux_match_errors(dd::IMAS.dd, par::FUSEparameters__ActorFluxMatcher, norms::Vector{Float64}, prog::Any)
    cp1d = dd.core_profiles.profiles_1d[]

    total_sources = IMAS.total_sources(dd.core_sources, cp1d; fields=[:total_ion_power_inside, :power_inside, :particles_inside, :torque_tor_inside])
    total_fluxes = IMAS.total_fluxes(dd.core_transport)

    cs_gridpoints = [argmin(abs.(rho_x .- total_sources.grid.rho_tor_norm)) for rho_x in par.rho_transport]
    cf_gridpoints = [argmin(abs.(rho_x .- total_fluxes.grid_flux.rho_tor_norm)) for rho_x in par.rho_transport]

    errors = Float64[]

    n = 0
    if par.evolve_Ti == :flux_match
        n += 1
        target = total_sources.total_ion_power_inside[cs_gridpoints] ./ total_sources.grid.surface[cs_gridpoints]
        output = total_fluxes.total_ion_energy.flux[cf_gridpoints]
        append!(errors, error_transformation(target, output, norms[n]))
    end

    if par.evolve_Te == :flux_match
        n += 1
        target = total_sources.electrons.power_inside[cs_gridpoints] ./ total_sources.grid.surface[cs_gridpoints]
        output = total_fluxes.electrons.energy.flux[cf_gridpoints]
        append!(errors, error_transformation(target, output, norms[n]))
    end

    if par.evolve_rotation == :flux_match
        n += 1
        target = total_sources.torque_tor_inside[cs_gridpoints] ./ total_sources.grid.surface[cs_gridpoints]
        output = total_fluxes.momentum_tor.flux[cf_gridpoints]
        append!(errors, error_transformation(target, output, norms[n]))
    end

    evolve_densities = evolve_densities_dictionary(cp1d, par)
    if !isempty(evolve_densities)
        if evolve_densities[:electrons] == :flux_match
            n += 1
            target = total_sources.electrons.particles_inside[cs_gridpoints] ./ total_sources.grid.surface[cs_gridpoints]
            output = total_fluxes.electrons.particles.flux[cf_gridpoints]
            append!(errors, error_transformation(target, output, norms[n]))
        end
        for ion in cp1d.ion
            if evolve_densities[Symbol(ion.label)] == :flux_match
                error("This is currently not working will fix later")
            end
        end
    end

    if prog !== nothing
        ProgressMeter.next!(prog; showvalues=progress_ActorFluxMatcher(dd, norm(errors)))
    end

    return errors
end

function progress_ActorFluxMatcher(dd::IMAS.dd, error::Float64)
    cp1d = dd.core_profiles.profiles_1d[]
    tmp = [
        ("         error", error),
        ("  Pfusion [MW]", IMAS.fusion_power(cp1d) / 1E6),
        ("     Ti0 [keV]", cp1d.ion[1].temperature[1] / 1E3),
        ("     Te0 [keV]", cp1d.electrons.temperature[1] / 1E3),
        ("ne0 [10²⁰ m⁻³]", cp1d.electrons.density_thermal[1] / 1E20)
    ]
    return tuple(tmp...)
end

function evolve_densities_dictionary(cp1d::IMAS.core_profiles__profiles_1d, par::FUSEparameters__ActorFluxMatcher)
    if par.evolve_densities == :fixed
        return setup_density_evolution_fixed(cp1d)
    elseif par.evolve_densities == :flux_match
        return setup_density_evolution_electron_flux_match_rest_ne_scale(cp1d)
    elseif typeof(par.evolve_densities) <: AbstractDict
        return par.evolve_densities
    else
        error(
            "act.ActorFluxMatcher.evolve_densities cannot be `$(repr(par.evolve_densities))`. Use `:fixed`, `:flux_match` or a dictionary specifiying [:quasi_neutrality, :match_ne_scale, :fixed, :flux_match] for each specie"
        )
    end
end

"""
    pack_z_profiles(cp1d::IMAS.core_profiles__profiles_1d, par::FUSEparameters__ActorFluxMatcher)

Packs the z_profiles based on evolution parameters

NOTE: the order for packing and unpacking is always: [Ti, Te, Rotation, ne, nis...]
"""
function pack_z_profiles(cp1d::IMAS.core_profiles__profiles_1d, par::FUSEparameters__ActorFluxMatcher)
    cp_gridpoints = [argmin(abs.(rho_x .- cp1d.grid.rho_tor_norm)) for rho_x in par.rho_transport]

    z_profiles = Float64[]

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

    evolve_densities = evolve_densities_dictionary(cp1d, par)
    if !isempty(evolve_densities)
        check_evolve_densities(cp1d, evolve_densities)
        if evolve_densities[:electrons] == :flux_match
            z_ne = IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal)[cp_gridpoints]
            append!(z_profiles, z_ne)
        end
        for ion in cp1d.ion
            if evolve_densities[Symbol(ion.label)] == :flux_match
                z_ni = IMAS.calcz(cp1d.grid.rho_tor_norm, ion.density_thermal)[cp_gridpoints]
                append!(z_profiles, z_ni)
            end
        end
    end

    return z_profiles
end

"""
    unpack_z_profiles(
        cp1d::IMAS.core_profiles__profiles_1d,
        par::FUSEparameters__ActorFluxMatcher,
        z_profiles::AbstractVector{<:Real})

Unpacks z_profiles based on evolution parameters

NOTE: The order for packing and unpacking is always: [Ti, Te, Rotation, ne, nis...]
"""
function unpack_z_profiles(
    cp1d::IMAS.core_profiles__profiles_1d,
    par::FUSEparameters__ActorFluxMatcher,
    z_profiles::AbstractVector{<:Real})

    cp_rho_transport = [cp1d.grid.rho_tor_norm[argmin(abs.(rho_x .- cp1d.grid.rho_tor_norm))] for rho_x in par.rho_transport]

    N = length(par.rho_transport)
    counter = 0

    if par.evolve_Ti == :flux_match
        Ti_new = IMAS.profile_from_z_transport(cp1d.ion[1].temperature, cp1d.grid.rho_tor_norm, cp_rho_transport, z_profiles[counter+1:counter+N])
        counter += N
        for ion in cp1d.ion
            ion.temperature = Ti_new
        end
    end

    if par.evolve_Te == :flux_match
        cp1d.electrons.temperature = IMAS.profile_from_z_transport(cp1d.electrons.temperature, cp1d.grid.rho_tor_norm, cp_rho_transport, z_profiles[counter+1:counter+N])
        counter += N
    end

    if par.evolve_rotation == :flux_match
        cp1d.rotation_frequency_tor_sonic =
            IMAS.profile_from_z_transport(cp1d.rotation_frequency_tor_sonic .+ 1, cp1d.grid.rho_tor_norm, cp_rho_transport, z_profiles[counter+1:counter+N])
        counter += N
    end

    evolve_densities = evolve_densities_dictionary(cp1d, par)
    if !isempty(evolve_densities)
        if evolve_densities[:electrons] == :flux_match
            cp1d.electrons.density_thermal =
                IMAS.profile_from_z_transport(cp1d.electrons.density_thermal, cp1d.grid.rho_tor_norm, cp_rho_transport, z_profiles[counter+1:counter+N])
            counter += N
        end
        z_ne = IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal)
        for ion in cp1d.ion
            if evolve_densities[Symbol(ion.label)] == :flux_match
                ion.density_thermal = IMAS.profile_from_z_transport(ion.density_thermal, cp1d.grid.rho_tor_norm, cp_rho_transport, z_profiles[counter+1:counter+N])
                counter += N
            elseif evolve_densities[Symbol(ion.label)] == :match_ne_scale
                ion.density_thermal = IMAS.profile_from_z_transport(ion.density_thermal, cp1d.grid.rho_tor_norm, cp1d.grid.rho_tor_norm, z_ne)
            end
        end
    end

    # Ensure Quasi neutrality if you are evolving densities
    if !IMAS.is_quasi_neutral(cp1d) && !isempty(evolve_densities) && any((evo != :fixed for evo in values(evolve_densities)))
        q_specie = Symbol[i for (i, evolve) in evolve_densities if evolve == :quasi_neutrality]
        @assert q_specie != Symbol[] "no quasi neutrality specie while quasi neutrality is broken: $(evolve_densities)"
        @assert length(q_specie) == 1 "only one specie can be set as :quasi_neutrality: $(evolve_densities)"
        IMAS.enforce_quasi_neutrality!(cp1d, q_specie[1])
    end

    return cp1d
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
    setup_density_evolution_electron_flux_match_rest_ne_scale(cp1d::IMAS.core_profiles__profiles_1d)

Sets up the evolve_density dict to evolve only ne and keep the rest matching the ne_scale lengths
"""
function setup_density_evolution_electron_flux_match_rest_ne_scale(cp1d::IMAS.core_profiles__profiles_1d)
    dd_thermal = [item[2] for item in IMAS.species(cp1d; only_electrons_ions=:ions, only_thermal_fast=:thermal)]
    dd_fast = [item[2] for item in IMAS.species(cp1d; only_electrons_ions=:all, only_thermal_fast=:fast)]
    quasi_neutrality_specie = :D
    if :DT ∈ dd_thermal
        quasi_neutrality_specie = :DT
    end
    return evolve_densities_dict_creation([:electrons], dd_fast, dd_thermal; quasi_neutrality_specie)
end

function setup_density_evolution_electron_flux_match_rest_ne_scale(dd::IMAS.dd)
    return setup_density_evolution_electron_flux_match_rest_ne_scale(dd.core_profiles.profiles_1d[])
end


"""
    setup_density_evolution_fixed(cp1d::IMAS.core_profiles__profiles_1d)

Sets up the evolve_density dict to keep all species fixed
"""
function setup_density_evolution_fixed(cp1d::IMAS.core_profiles__profiles_1d)
    s = [item[2] for item in IMAS.species(cp1d)]
    return evolve_densities_dict_creation(Symbol[], s, Symbol[]; quasi_neutrality_specie=false)
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
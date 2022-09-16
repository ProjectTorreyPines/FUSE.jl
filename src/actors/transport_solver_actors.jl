import NLsolve: nlsolve

#= ========================== =#
#     transport solver actor   #
#= ========================== =#

mutable struct ActorTransportSolver <: PlasmaAbstractActor
    dd::IMAS.dd
    par::ParametersActor
end

function ParametersActor(::Type{Val{:ActorTransportSolver}})
    par = ParametersActor(nothing)
    par.evolve_Ti = Switch([:flux_match, :fixed], "", "How to evolve the average ion temperature "; default=:flux_match)
    par.evolve_Te = Switch([:flux_match, :fixed], "", "How to evolve the electron temperature"; default=:flux_match)
    par.evolve_densities = Entry(Union{Dict,Symbol}, "", "OrderedDict to specify which ion species are evolved, kept constant or used for quasi neutarlity"; default=:fixed)
    par.evolve_rotation = Switch([:flux_match, :fixed], "", "How to evolve the electron temperature"; default=:fixed)
    par.rho_transport = Entry(Vector{Float64}, "", "Rho transport grid")
    par.do_plot = Entry(Bool, "", "plots the flux matching"; default=false)
    return par
end

"""
    ActorTransportSolver(dd::IMAS.dd, act::ParametersAllActors; kw...)


"""
function ActorTransportSolver(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorTransportSolver(kw...)
    actor = ActorTransportSolver(dd, par)
    step(actor, act)
    finalize(actor)
    return actor
end

"""
function ActorTransportSolver(dd::IMAS.dd, par::ParametersActor; kw...)
    par = par(kw...)
    return ActorTransportSolver(dd, par)
end
"""

function step(actor::ActorTransportSolver, act::ParametersAllActors; iterations::Int=50, factor=1.0)
    dd = actor.dd
    par = actor.par

    z_init = pack_z_profiles(dd, par)
    dd_before = deepcopy(dd)
    @show length(z_init) == length(get_flux_match_error(dd, act, z_init))
    res = nlsolve(z -> get_flux_match_error(dd, act, z), z_init, show_trace=true, factor=factor)#, iterations=iterations, factor=factor)
#    dd = unpack_z_profiles(dd_before, par, res.minimum)
    dd = unpack_z_profiles(dd_before, par, res.minimum)
#    get_flux_match_error(dd,act, res.minimum)
    get_flux_match_error(dd,act, res.minimum)

    actor.dd = dd

    if par.do_plot
        display(plot(dd.core_transport))
        plot(dd_before.core_profiles,label="before")
        display(plot!(dd.core_profiles,label="after"))
    end

    return actor
end


function get_flux_match_error(dd::IMAS.dd, act::ParametersAllActors, z_profiles::AbstractVector{<:Float64})
    par = act.ActorTransportSolver
    unpack_z_profiles(dd, par, z_profiles) # modify dd with new z_profiles
    act.ActorTGLF.tglf_model = :tglfnn
    act.ActorTGLF.rho_transport = par.rho_transport
    act.ActorNeoclassical.rho_transport = par.rho_transport
    FUSE.ActorTGLF(dd, act)
    FUSE.ActorNeoclassical(dd, act)
    IMAS.sources!(dd)
    return pack_flux_match_errors(dd, par)
end


function pack_flux_match_errors(dd, par)
    total_sources = IMAS.total_sources(dd)
    total_fluxes = IMAS.total_fluxes(dd)

    cp1d = dd.core_profiles.profiles_1d[]
    cs_gridpoints = [argmin((rho_x .- total_sources.grid.rho_tor_norm).^2) for rho_x in par.rho_transport]
    ct_gridpoints = [argmin((rho_x .- total_fluxes.grid_flux.rho_tor_norm).^2) for rho_x in par.rho_transport]
    surface = total_sources.grid.surface[cs_gridpoints]

    energy_norm = 1e4
    particles_norm = 1e19
    momentum_norm = 0.1

    error = Vector{Float64}()

    if par.evolve_Ti == :flux_match
        append!(error, (total_sources.total_ion_power_inside[cs_gridpoints] ./ surface .- total_fluxes.total_ion_energy.flux[ct_gridpoints]) ./ energy_norm)
    end

    if par.evolve_Te == :flux_match
        append!(error, (total_sources.electrons.power_inside[cs_gridpoints] ./ surface .- total_fluxes.electrons.energy.flux[ct_gridpoints]) ./ energy_norm)
    end

    if par.evolve_rotation == :flux_match
        append!(error, (total_sources.torque_tor_inside[cs_gridpoints] ./ surface .- total_fluxes.momentum_tor.flux[ct_gridpoints]) ./ momentum_norm)
    end

    if par.evolve_densities != :fixed
        if par.evolve_densities[:electrons] == :flux_match
            append!(error, (total_sources.electrons.particles_inside[cs_gridpoints] ./ surface .- total_fluxes.electrons.particles.flux[ct_gridpoints]) ./ particles_norm)
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
Note the order for packing and unpacking is always:
[Ti, Te, Rotation, ne, nis...]
"""

function pack_z_profiles(dd, par)
    cp1d = dd.core_profiles.profiles_1d[]
    z_profiles = Vector{Float64}()
    cp_gridpoints = [argmin((rho_x .- cp1d.grid.rho_tor_norm).^2) for rho_x in par.rho_transport]

    if par.evolve_Ti == :flux_match
        append!(z_profiles, (IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.ion[1].temperature)[cp_gridpoints]))
    end

    if par.evolve_Te == :flux_match
        append!(z_profiles, (IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.electrons.temperature)[cp_gridpoints]))
    end

    if par.evolve_rotation == :flux_match
        append!(z_profiles, (IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.rotation_frequency_tor_sonic)[cp_gridpoints]))
    end

    if par.evolve_densities != :fixed
        if par.evolve_densities[:electrons] == :flux_match
            append!(z_profiles, (IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.electrons.density_thermal)[cp_gridpoints]))
        end
        for ion in cp1d.ion
            if par.evolve_densities[Symbol(ion.label)] == :flux_match
                append!(z_profiles, (IMAS.calcz(cp1d.grid.rho_tor_norm, ion.density_thermal)[cp_gridpoints]))
            end
        end
    end
    return z_profiles
end

function unpack_z_profiles(dd, par, z_profiles)
    cp1d = dd.core_profiles.profiles_1d[]
    rho_transport = par.rho_transport
    counter = 0
    N = Int(length(rho_transport))

    if par.evolve_Ti == :flux_match
        Ti_new = IMAS.profile_from_z_transport(cp1d.ion[1].temperature, cp1d.grid.rho_tor_norm, rho_transport, z_profiles[counter+1:counter+N])
        counter += N
        for ion in cp1d.ion
            ion.temperature = Ti_new
        end
    end

    if par.evolve_Te == :flux_match
        cp1d.electrons.temperature = IMAS.profile_from_z_transport(cp1d.electrons.temperature, cp1d.grid.rho_tor_norm, rho_transport, z_profiles[counter+1:counter+N])
        counter += 4
    end

    if par.evolve_rotation == :flux_match
        cp1d.rotation_frequency_tor_sonic = IMAS.profile_from_z_transport(cp1d.rotation_frequency_tor_sonic, cp1d.grid.rho_tor_norm, rho_transport, z_profiles[counter+1:counter+N])
        counter += 4
    end

    if par.evolve_densities != :fixed
        if par.evolve_densities[:electrons] == :flux_match
            cp1d.electrons.density_thermal = IMAS.profile_from_z_transport(cp1d.electrons.density_thermal, cp1d.grid.rho_tor_norm, rho_transport, z_profiles[counter+1:counter+N])
            counter += 4
        end
        for ion in cp1d.ion
            if par.evolve_densities[Symbol(ion.label)] == :flux_match
                ion.density_thermal = IMAS.profile_from_z_transport(ion.density_thermal, cp1d.grid.rho_tor_norm, rho_transport, z_profiles[counter+1:counter+N])
                counter += 4
            end
        end
    end
    return dd
end

function check_evolve_densities(dd::IMAS.dd, evolve_densities::Dict)
    dd_species = vcat([Symbol(ion.label) for ion in dd.core_profiles.profiles_1d[].ion], :electrons)
    # Check if evolve_densities contains all of dd species
    @assert sort([i for (i, evolve) in evolve_densities]) == sort(dd_species) "Not all species are accounted for in the evolve_densities dict : $(sort([i for (i,j) in evolve_densities])) , dd_species : $(sort(dd_species)) ,"
    # Check if there is 1 quasi_neutality specie
    @assert length([i for (i, evolve) in evolve_densities if evolve == :quasi_neutality]) < 2 "Only one specie can be used for quasi neutality matching"
    # Check if there is a source for the flux_match specie(s)

    return true
end


"""


Sets up the evolve_density dict to evolve only ne and keep the rest matching the ne_scale lengths
"""
function setup_density_evolution(dd)
    dd_species = [Symbol(ion.label) for ion in dd.core_profiles.profiles_1d[].ion]

    return evolve_densities_dict_creation([:electrons], [], dd_species)
end

function evolve_densities_dict_creation(flux_match_species::Vector, fixed_species::Vector, match_ne_scale_species::Vector; quasi_neutality_specie::Union{Symbol,Bool}=false)
    parse_list = vcat([[sym, :flux_match] for sym in flux_match_species], [[sym, :match_ne_scale] for sym in match_ne_scale_species], [[sym, :fixed] for sym in fixed_species])
    if isa(quasi_neutality_specie, Symbol)
        parse_list = vcat(parse_list, [[quasi_neutality_specie, :quasi_neutality]])
    end

    return Dict(sym => evolve for (sym, evolve) in parse_list)
end

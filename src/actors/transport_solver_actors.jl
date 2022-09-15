#using NLsolve

#= ========================== =#
#     transport solver actor   #
#= ========================== =#

mutable struct ActorTransportSolver <: PlasmaAbstractActor
    dd::IMAS.dd
    par::ParametersActor
end

function ParametersActor(::Type{Val{:ActorTransportSolver}})
    par = ParametersActor(nothing)
    par.evolve_Ti = Switch([:flux_match, :fixed], "", "How to evolve the average ion temperature "; default=:evolve)
    par.evolve_Te = Switch([:flux_match, :fixed], "", "How to evolve the electron temperature"; default=:evolve)
    par.evolve_densities = Entry(Union{Dict,Symbol}, "", "OrderedDict to specify which ion species are evolved, kept constant or used for quasi neutarlity"; default=:fixed)
    par.evolve_rotation = Switch([:flux_match, :fixed], "", "How to evolve the electron temperature"; default=:fixed)
    par.rho_transport = Entry(Vector{Real}, "", "Rho transport grid")
    return par
end

"""
    ActorTransportSolver(dd::IMAS.dd, act::ParametersAllActors; kw...)


"""
function ActorTransportSolver(dd::IMAS.dd, act::ParametersAllActors; kw...)
    par = act.ActorTransportSolver(kw...)
    actor = ActorTransportSolver(dd, par)
    step(actor)
    finalize(actor)
    return actor
end

function ActorTransportSolver(dd::IMAS.dd, par::ParametersActor; kw...)
    par = par(kw...)
    return ActorTransportSolver(dd, par)
end

function step(actor::ActorTransportSolver; iterations::Int=50, factor=1.0)
    dd = actor.dd
    par = actor.par

    z_init = pack_z_profiles(dd, par)
    dd_before = deepcopy(dd)

    res = nlsolve(z -> get_flux_match_error(dd, par.rho_transport, z), z_init, show_trace=true, iterations=iterations, factor=factor)
    dd = unpack(dd_before, res.minimum, rho_transport)
    get_flux_match_error(dd, rho_transport, res.minimum)
    actor.dd = dd
    return actor
end


function get_flux_match_error(dd::IMAS.dd, par::ParametersActor, z_profiles::AbstractVector{<:Real})
    unpack(dd, z_profiles, rho_transport) # modify dd with new z_profiles
    act.ActorTGLF.tglf_model = :tglfnn
    act.ActorTGLF.rho_tglf = par.rho_transport
    FUSE.ActorTGLF(dd, act)

    IMAS.sources!(dd)

    pack_flux_match_errors(dd, par)

end


function pack_flux_match_errors(dd, par)
    total_sources = IMAS.total_sources(dd)
    total_fluxes = IMAS.total_fluxes(dd)

    energy_norm = 1e4
    particles_norm = 1e19
    return error
end



"""
Note the order for packing and unpacking is always:
[Ti, Te, Rotation, ne, nis...]
"""

function pack_z_profiles(dd, par)
    cp1d = dd.core_profiles.profiles_1d[]
    z_profiles = Vector{Real}()

    if par.evolve_Ti
        append!(z_profiles, IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.ion[1].temperature))
    end

    if par.evolve_Te
        append!(z_profiles, IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.electrons.temperature))
    end

    if par.evolve_rotation
        append!(z_profiles, IMAS.calc_z(cp1d.grid.rho_tor_norm, cp1d.rotation_frequency_tor_sonic))
    end

    if par.evolve_densities != :fixed
        if par.evolve_densities[:electrons] == :flux_match
            append!(z_profiles, cp1d.electrons.density_thermal)
        end
        for ion in cp1d.ion
            if par.evolve_densities[Symbol(ion.label)] == :flux_match
                append!(z_profiles, ion.density_thermal)
            end
        end
    end
    return z_profiles
end

function unpack(dd, z_profiles, rho_transport)
    cp1d = dd.core_profiles.profiles_1d[]
    counter = 0
    N = length(z_profiles) / length(rho_transport)

    if par.evolve_Ti
        Ti_new = IMAS.profile_from_z_transport(cp1d.ion[1].temperature, cp1d.grid.rho_tor_norm, rho_transport, z_profiles[counter+1:counter+N])
        counter += N
        for ion in cp1d.ion
            ion.temperature = Ti_new
        end
    end

    if par.evolve_Te
        cp1d.electrons.temperature = IMAS.profile_from_z_transport(cp1d.electrons.temperature, cp1d.grid.rho_tor_norm, rho_transport, z_profiles[counter+1:counter+N])
        counter += 4
    end

    if par.evolve_rotation
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

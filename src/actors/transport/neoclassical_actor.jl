import NEO
import GACODE
import Interpolations
import MillerExtendedHarmonic

#= ================= =#
#  ActorNeoclassical  #
#= ================= =#
Base.@kwdef mutable struct FUSEparameters__ActorNeoclassical{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    model::Switch{Symbol} = Switch{Symbol}([:changhinton, :neo, :hirshmansigmar_facit], "-", "Neoclassical model to run"; default=:hirshmansigmar_facit)
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "rho_tor_norm values to compute neoclassical fluxes on"; default=0.25:0.1:0.85)
end

mutable struct ActorNeoclassical{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorNeoclassical{P}}
    input_neos::Vector{<:NEO.InputNEO}
    flux_solutions::Vector{<:GACODE.FluxSolution}
    equilibrium_geometry::Union{NEO.EquilibriumGeometry,Missing}
    facit_output::Union{Dict{String, NEO.FACIToutput},Missing}
end

"""
    ActorNeoclassical(dd::IMAS.dd, act::ParametersAllActors; kw...)

Evaluates the neoclassical transport fluxes
"""
function ActorNeoclassical(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorNeoclassical(dd, act.ActorNeoclassical; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function ActorNeoclassical(dd::IMAS.dd, par::FUSEparameters__ActorNeoclassical; kw...)
    logging_actor_init(ActorNeoclassical)
    par = OverrideParameters(par; kw...)
    return ActorNeoclassical(dd, par, NEO.InputNEO[], GACODE.FluxSolution[], missing, missing)
end

"""
    _step(actor::ActorNeoclassical)

Runs ActorNeoclassical to evaluate the neoclassical transport flux on a vector of gridpoints
"""
function _step(actor::ActorNeoclassical)
    par = actor.par
    dd = actor.dd

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]
    rho_cp = cp1d.grid.rho_tor_norm

    if par.model == :changhinton
        actor.flux_solutions = [NEO.changhinton(eqt, cp1d, rho, 1) for rho in par.rho_transport]

    elseif par.model == :neo
        gridpoint_cps = [argmin_abs(rho_cp, rho) for rho in par.rho_transport]
        actor.input_neos = [NEO.InputNEO(eqt, cp1d, i) for (idx, i) in enumerate(gridpoint_cps)]
        actor.flux_solutions = asyncmap(input_neo -> NEO.run_neo(input_neo), actor.input_neos)

    elseif par.model == :hirshmansigmar_facit
        gridpoint_cps = [argmin_abs(rho_cp, rho) for rho in par.rho_transport]
        if ismissing(actor.equilibrium_geometry) || actor.equilibrium_geometry.time != eqt.time
            actor.equilibrium_geometry = NEO.get_equilibrium_geometry(eqt, cp1d)#, gridpoint_cps)
        end
        parameter_matrices = NEO.get_plasma_profiles(eqt, cp1d)
        rho_s = GACODE.rho_s(cp1d, eqt)
        rmin = GACODE.r_min_core_profiles(eqt.profiles_1d, cp1d.grid.rho_tor_norm)
        actor.flux_solutions = map(gridpoint_cp -> NEO.hirshmansigmar(gridpoint_cp, eqt, cp1d, parameter_matrices, actor.equilibrium_geometry; rho_s, rmin), gridpoint_cps)
        
        facit_full_geometry = true 
        # impurity fluxes with FACIT
        outputs_dict = Dict{String, NEO.FACIToutput}()

        for i in 1:length(cp1d.ion)-1
            for j in i+1:length(cp1d.ion)
                species1 = cp1d.ion[i]
                species2 = cp1d.ion[j]

                z = species2.element[1].z_n
                if z ≤ 18
                    facit_rotation_model = 0
                else
                    facit_rotation_model = 2
                end

                key = species1.label * "+" * species2.label
                facit_input = prepare_facit(dd, facit_rotation_model, facit_full_geometry, species1, species2)
                output = NEO.compute_transport(facit_input)

                outputs_dict[key] = output
            end
        end

        actor.facit_output = outputs_dict
    end

    return actor
end

"""
    _finalize(actor::ActorNeoclassical)

Writes ActorNeoclassical results to dd.core_transport
"""
function _finalize(actor::ActorNeoclassical)
    par = actor.par
    dd = actor.dd

    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]

    model = resize!(dd.core_transport.model, :neoclassical; wipe=false)
    m1d = resize!(model.profiles_1d)
    m1d.grid_flux.rho_tor_norm = par.rho_transport

    if par.model == :changhinton
        model.identifier.name = "Chang-Hinton"
        GACODE.flux_gacode_to_imas((:ion_energy_flux,), actor.flux_solutions, m1d, eqt, cp1d)

    elseif par.model == :neo
        model.identifier.name = "NEO"
        GACODE.flux_gacode_to_imas((:electron_energy_flux, :ion_energy_flux, :electron_particle_flux, :ion_particle_flux, :momentum_flux), actor.flux_solutions, m1d, eqt, cp1d)

    elseif par.model == :hirshmansigmar_facit
        model.identifier.name = "Hirshman-Sigmar + FACIT"
        GACODE.flux_gacode_to_imas((:electron_energy_flux, :ion_energy_flux, :electron_particle_flux, :ion_particle_flux), actor.flux_solutions, m1d, eqt, cp1d)

        hs_index = findfirst(model -> model.identifier.name == "Hirshman-Sigmar + FACIT", dd.core_transport.model)
  
        flux_by_species = Dict{String, Vector{Float64}}()

        for (key, result) in actor.facit_output
            first_label, second_label = split(key, "+")
            if haskey(flux_by_species, second_label)
                flux_by_species[second_label] .+= result.Flux_z
            else
                flux_by_species[second_label] = copy(result.Flux_z)
            end
        end

        rho = first(values(actor.facit_output)).rho
        for (label, total_flux_z) in flux_by_species
            interpolated_flux = IMAS.interp1d(rho, total_flux_z).(par.rho_transport)

            for ion in dd.core_transport.model[hs_index].profiles_1d[].ion
                if ion.label == label
                    ion.particles.flux = interpolated_flux
                end
            end
        end
    end

    return actor
end

function prepare_facit(dd::IMAS.dd, facit_rotation_model::Int, facit_full_geometry::Bool, species1::IMAS.core_profiles__profiles_1d___ion{Float64}, species2::IMAS.core_profiles__profiles_1d___ion{Float64})
    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]
    rmin = GACODE.r_min_core_profiles(eqt.profiles_1d, cp1d.grid.rho_tor_norm) ./ 1e2 
    a = rmin[end]

    function facit_interpolate(profile::Vector{Float64}, original_grid::Vector{Float64}, rmin_over_a::Vector{Float64})
        return IMAS.interp1d(original_grid, profile).(rmin_over_a)
    end

    rho = rmin ./ a
    theta = collect(range(0, 2π, length = length(rho)))
    Zimp = species2.z_ion
    Aimp = species2.element[1].a 
    Zi = species1.z_ion 
    Ai = species1.element[1].a 

    Ti = facit_interpolate(species1.temperature, cp1d.grid.rho_tor_norm, rho)
    Te = facit_interpolate(cp1d.electrons.temperature, cp1d.grid.rho_tor_norm, rho)
    Te_Ti = Te ./ Ti
    Ni = facit_interpolate(species1.density, cp1d.grid.rho_tor_norm, rho)
    Nimp = facit_interpolate(species2.density, cp1d.grid.rho_tor_norm, rho)

    if facit_rotation_model == 0
        Machi = zeros(length(rho))
    else 
        Vphi = cp1d.rotation_frequency_tor_sonic
        vth_i = sqrt.(2 .* Ti .* IMAS.mks.e ./ (Ai * IMAS.mks.m_p))
        Machi = Vphi ./ vth_i
    end

    Zeff = facit_interpolate(cp1d.zeff, cp1d.grid.rho_tor_norm, rho)
    gradTi = IMAS.gradient(rho .* a, Ti)
    gradNi = IMAS.gradient(rho .* a, Ni)
    gradNimp = IMAS.gradient(rho .* a, Nimp)

    invaspct = eqt.boundary.minor_radius / eqt.boundary.geometric_axis.r
    B0 = abs(eqt.global_quantities.vacuum_toroidal_field.b0)
    R0 = eqt.boundary.geometric_axis.r
    qmag = abs.(facit_interpolate(eqt.profiles_1d.q, eqt.profiles_1d.rho_tor_norm, rho))
    psi = facit_interpolate(eqt.profiles_1d.psi, eqt.profiles_1d.rho_tor_norm, rho)
    FV = facit_interpolate(eqt.profiles_1d.f, eqt.profiles_1d.psi, psi)
    dpsidx = IMAS.gradient(rho, psi)
    R, Z = FUSE.prepare_R_Z(dd, rho, theta)
    RV = R
    ZV = Z

    fj0 = NEO.FACITinput(rho, Zimp, Aimp, Zi, Ai, Ti, Ni, Nimp, Machi, Zeff, gradTi, gradNi, gradNimp, invaspct, B0, R0, qmag; 
        fsaout = false, rotation_model = facit_rotation_model, Te_Ti = Te_Ti, full_geom = facit_full_geometry, theta = theta, RV = RV, FV = FV, ZV = ZV, BV = missing, JV = missing, dpsidx = dpsidx, nat_asym = true)
    
    return fj0
end

function prepare_R_Z(dd::IMAS.dd, rho::Vector{Float64}, theta::Vector{Float64})
    eqt = dd.equilibrium.time_slice[]
    fw = IMAS.first_wall(dd.wall)
    surfaces = IMAS.trace_surfaces(eqt, fw.r, fw.z)

    function rho_geom_to_psi(surfaces::Vector{IMAS.FluxSurface})
        rmin_values = [fs.min_r for fs in surfaces]
        rmin_axis = minimum(rmin_values)
        rmin_edge = maximum(rmin_values)
        a = rmin_edge - rmin_axis
        
        psi_vals = [fs.psi for fs in surfaces]
        rho_geom = [(r - rmin_axis) / a for r in rmin_values]

        sorted = sortperm(rho_geom)
        rho_sorted = rho_geom[sorted]
        
        return Interpolations.LinearInterpolation(rho_sorted, psi_vals; extrapolation_bc=Interpolations.Line())
    end

    interp_psi = rho_geom_to_psi(surfaces)
    target_psi = interp_psi.(rho)

    sortf = sortperm(eqt.profiles_1d.f)
    f_sorted = eqt.profiles_1d.f[sortf]

    sortpsi = sortperm(eqt.profiles_1d.psi)
    psi_sorted = eqt.profiles_1d.psi[sortpsi]

    interp_f = Interpolations.LinearInterpolation(psi_sorted, f_sorted; extrapolation_bc=Interpolations.Line())

    target_f = interp_f.(target_psi)

    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    r, z, PSI_interpolant = IMAS.ψ_interpolant(eqt2d)
    RA = eqt.global_quantities.magnetic_axis.r
    ZA = eqt.global_quantities.magnetic_axis.z
    Br, Bz = IMAS.Br_Bz(eqt2d)

    selected_surfaces = IMAS.trace_surfaces(target_psi, target_f, r, z, eqt2d.psi, Br, Bz, PSI_interpolant, RA, ZA, fw.r, fw.z)

    function RZ_at_rtheta(coeffs::MillerExtendedHarmonic.MXH, theta::Vector{Float64})
        n = length(theta)
        R_rth = Vector{Float64}(undef, n)
        Z_rth = Vector{Float64}(undef, n)
    
        for i in eachindex(theta)
            R_rth[i] = MillerExtendedHarmonic.R_MXH(theta[i], coeffs)
            Z_rth[i] = MillerExtendedHarmonic.Z_MXH(theta[i], coeffs)
        end
    
        return R_rth, Z_rth
    end
    
    function miller_to_RZrtheta(surfaces::Vector{IMAS.FluxSurface}, theta::Vector{Float64})
        nthetas = length(theta)
        RV = zeros(length(surfaces), nthetas)
        ZV = zeros(length(surfaces), nthetas) 
        
        for (i, surface) in enumerate(surfaces)
            MillerExtendedHarmonic.reorder_flux_surface!(surface.r, surface.z)
            coeffs = MillerExtendedHarmonic.MXH(surface.r, surface.z)
    
            R, Z = RZ_at_rtheta(coeffs, theta)
    
            RV[i, :] .= R
            ZV[i, :] .= Z
        end
    
        return RV, ZV
    end

    RV, ZV = miller_to_RZrtheta(selected_surfaces, theta)

end
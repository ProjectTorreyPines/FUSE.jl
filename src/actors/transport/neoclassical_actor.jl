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
    model_bulk::Switch{Symbol} = Switch{Symbol}([:changhinton, :neo, :hirshmansigmar], "-", "Neoclassical model to run for bulk"; default=:hirshmansigmar)
    model_impurities::Switch{Symbol} = Switch{Symbol}([:facit, :hirshmansigmar], "-", "Neoclassical model to run for impurities"; default=:facit)
    facit_full_geometry::Entry{Bool} = Entry{Bool}("-", "Use full geometry (true) or circular approximation (false)"; default=true)
    facit_rotation_model::Entry{Int} = Entry{Int}("-", "Rotation model to use in FACIT"; default=0, check=x -> @assert 0 ≤ x ≤ 2 "rotation_model options are 0, 1, or 2")
    rho_transport::Entry{AbstractVector{T}} = Entry{AbstractVector{T}}("-", "rho_tor_norm values to compute neoclassical fluxes on"; default=0.25:0.1:0.85)
end

mutable struct ActorNeoclassical{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorNeoclassical{P}}
    input_neos::Vector{<:NEO.InputNEO}
    flux_solutions::Vector{<:GACODE.FluxSolution}
    equilibrium_geometry::Union{NEO.EquilibriumGeometry,Missing}
    facit_output::Union{NEO.FACIToutput,Missing}
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

    if par.model_bulk == :changhinton
        actor.flux_solutions = [NEO.changhinton(eqt, cp1d, rho, 1) for rho in par.rho_transport]

    elseif par.model_bulk == :neo
        gridpoint_cps = [argmin_abs(rho_cp, rho) for rho in par.rho_transport]
        actor.input_neos = [NEO.InputNEO(eqt, cp1d, i) for (idx, i) in enumerate(gridpoint_cps)]
        actor.flux_solutions = asyncmap(input_neo -> NEO.run_neo(input_neo), actor.input_neos)

    elseif par.model_bulk == :hirshmansigmar
        gridpoint_cps = [argmin_abs(rho_cp, rho) for rho in par.rho_transport]
        if ismissing(actor.equilibrium_geometry) || actor.equilibrium_geometry.time != eqt.time
            actor.equilibrium_geometry = NEO.get_equilibrium_geometry(eqt, cp1d)#, gridpoint_cps)
        end
        parameter_matrices = NEO.get_plasma_profiles(eqt, cp1d)
        actor.flux_solutions = map(gridpoint_cp -> NEO.hirshmansigmar(gridpoint_cp, eqt, cp1d, parameter_matrices, actor.equilibrium_geometry), gridpoint_cps)
        
        if par.model_impurities == :facit 
            fj0 = prepare_facit(actor)
            actor.facit_output = NEO.compute_transport(fj0)
        end 
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

    if par.model_bulk == :changhinton
        model.identifier.name = "Chang-Hinton"
        GACODE.flux_gacode_to_imas((:ion_energy_flux,), actor.flux_solutions, m1d, eqt, cp1d)

    elseif par.model_bulk == :neo
        model.identifier.name = "NEO"
        GACODE.flux_gacode_to_imas((:electron_energy_flux, :ion_energy_flux, :electron_particle_flux, :ion_particle_flux, :momentum_flux), actor.flux_solutions, m1d, eqt, cp1d)

    elseif par.model_bulk == :hirshmansigmar
        model.identifier.name = "Hirshman-Sigmar"
        GACODE.flux_gacode_to_imas((:electron_energy_flux, :ion_energy_flux, :electron_particle_flux, :ion_particle_flux), actor.flux_solutions, m1d, eqt, cp1d)

        if par.model_impurities == :facit
            # interpolate onto the same grid as hirshman sigmar 
            flux_rho_transport = IMAS.interp1d(actor.facit_output.rho, actor.facit_output.Flux_z).(par.rho_transport)
            hs_index = findfirst(model -> model.identifier.name == "Hirshman-Sigmar", dd.core_transport.model)
            dd.core_transport.model[hs_index].identifier.name = "Hirshman-Sigmar + FACIT"
            dd.core_transport.model[hs_index].profiles_1d[].ion[end].particles.flux = flux_rho_transport
        end
    end

    return actor
end

function prepare_facit(actor::ActorNeoclassical)
    dd = actor.dd
    par = actor.par

    cp1d = dd.core_profiles.profiles_1d[]
    eqt = dd.equilibrium.time_slice[]
    rmin = GACODE.r_min_core_profiles(eqt.profiles_1d, cp1d.grid.rho_tor_norm) ./ 1e2 
    a = rmin[end]

    function facit_interpolate(profile::Vector{Float64}, original_grid::Vector{Float64}, rmin_over_a::Vector{Float64})
        return IMAS.interp1d(original_grid, profile).(rmin_over_a)
    end

    rho = rmin ./ a
    theta = range(0, 2π, length = length(rho))
    Zimp = cp1d.ion[end].z_ion
    Aimp = cp1d.ion[end].element[1].a
    Zi = cp1d.ion[1].z_ion
    Ai = cp1d.ion[1].element[1].a
    Ti = facit_interpolate(cp1d.ion[1].temperature, cp1d.grid.rho_tor_norm, rho)
    ne = facit_interpolate(cp1d.electrons.density, cp1d.grid.rho_tor_norm, rho)
    Ni = facit_interpolate(cp1d.ion[1].density, cp1d.grid.rho_tor_norm, rho)
    Nimp = facit_interpolate(cp1d.ion[end].density, cp1d.grid.rho_tor_norm, rho)
    Mi_core = 0.35
    Mi_edge = 0.05

    Machi = (Mi_core - Mi_edge) * (1 .- rho.^2) .+ Mi_edge
    Zeff = facit_interpolate(cp1d.zeff, cp1d.grid.rho_tor_norm, rho)
    gradTi = -IMAS.calc_z(rho, cp1d.ion[1].temperature, :backward)
    gradNi = -IMAS.calc_z(rho, cp1d.ion[1].density, :backward)
    gradNimp = -IMAS.calc_z(rho, cp1d.ion[end].density, :backward)
    invaspct = eqt.boundary.minor_radius / eqt.boundary.geometric_axis.r
    B0 = -eqt.global_quantities.vacuum_toroidal_field.b0
    R0 = eqt.boundary.geometric_axis.r
    qmag = -facit_interpolate(eqt.profiles_1d.q, eqt.profiles_1d.rho_tor_norm, rho)
    FV = facit_interpolate(eqt.profiles_1d.f, eqt.profiles_1d.rho_tor_norm, rho)
    psi = facit_interpolate(eqt.profiles_1d.psi, eqt.profiles_1d.rho_tor, rho)
    dpsidx = IMAS.gradient(rho, psi)

    R, Z = FUSE.prepare_R_Z(dd, rho, collect(theta))
    RV = R
    ZV = Z

    fj0 = NEO.FACITinput(rho, Zimp, Aimp, Zi, Ai, Ti, Ni, Nimp, Machi, Zeff, gradTi, gradNi, gradNimp, invaspct, B0, R0, qmag; 
        fsaout = true, rotation_model = par.facit_rotation_model, full_geom = par.facit_full_geometry, RV = RV, FV = FV, ZV = ZV, BV = missing, JV = missing, dpsidx = dpsidx, nat_asym = true)
    
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

    interp_f = Interpolations.LinearInterpolation(eqt.profiles_1d.psi, eqt.profiles_1d.f; extrapolation_bc=Interpolations.Line())
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
    
            RV[:, i] .= R
            ZV[i, :] .= Z
        end
    
        return RV, ZV
    end

    RV, ZV = miller_to_RZrtheta(selected_surfaces, theta)

end
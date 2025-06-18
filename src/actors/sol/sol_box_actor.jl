#= =========== =#
#  ActorSOLBox  #
#= =========== =#

Base.@kwdef mutable struct FUSEparameters__ActorSOLBox{T<:Real} <: ParametersActor{T}
    _parent::WeakRef = WeakRef(nothing)
    _name::Symbol = :not_set
    _time::Float64 = NaN
    Te_t::Entry{T} = Entry{T}("eV", "Input electron temperature at the target"; default=10.0)
    Ti_t::Entry{T} = Entry{T}("eV", "Input ion temperature at the target"; default=10.0)
    frac_cond::Entry{T} = Entry{T}("-", "Fraction of power carried by electron conduction"; default=0.7)
    frac_mom::Entry{T} = Entry{T}("-", "Fraction of momentum lost due to collisions with neutrals, atomic processes and viscous forces"; default=0.5)
    κ0_e::Entry{T} = Entry{T}("-", "Coefficient of electron conductivity"; default=2000.0)
    κ0_i::Entry{T} = Entry{T}("-", "Coefficient of ion conductivity"; default=60.0)
    recycling_coeff_i::Entry{T} = Entry{T}("-", "Ion particle recycling coefficient"; default=0.99)
    recycling_coeff_e::Entry{T} = Entry{T}("-", "Electron particle recycling coefficient"; default=0.99)
    λq::Entry{T} = Entry{T}("m", "Width of the flux tube"; default=Inf)
    verbose::Entry{Bool} = act_common_parameters(; verbose=false)
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

Base.@kwdef mutable struct SOLBox{T<:Real}
    single_null::Bool = false
    κ0_i::T = T(NaN) # Coefficient of ion conductivity
    recycling_coeff_i::T = T(NaN) # Ion particle recycling coefficient
    Ti_t::T = T(NaN) # Input ion temperature at the target [eV]
    mass_ion::T = T(NaN) # Ion mass [AMU]
    κ0_e::T = T(NaN) # Coefficient of electron conductivity
    recycling_coeff_e::T = T(NaN) # Electron particle recycling coefficient
    Te_t::T = T(NaN) # Input electron temperature at the target [eV]
    frac_cond::T = T(NaN) # Fraction of power carried by electron conduction
    frac_mom::T = T(NaN) # Fraction of momentum lost due to collisions with neutrals, atomic processes and viscous forces
    λq::T = T(NaN) # Width of the flux tube
    qpar_i::T = T(NaN) # Upstream parallel ion heat flux [W m^-2]
    qpar_e::T = T(NaN) # Upstream parallel electron heat flux [W m^-2]
    Γ_i::T = T(NaN) # Upstream ion particle flux [ptcles m^-2 s^-1]
    Γ_e::T = T(NaN)  # Upstream electron particle flux [ptcles m^-2 s^-1]
    Te_u::T = T(NaN) # Electron temperature at the upstream (separatrix) [eV]
    Ti_u::T = T(NaN) # Ion temperature at the upstream (separatrix) [eV]
    ne_t::T = T(NaN) # Electron density at the target [pcles m^-2 s^-1]
    ni_t::T = T(NaN) # Ion density at the target [pcles m^-2 s^-1]
    ne_u::T = T(NaN) # Electron density at the upstream (separatrix) [pcles m^-2 s^-1]
    ni_u::T = T(NaN) # Ion density at the upstream (separatrix) [pcles m^-2 s^-1]
    cst::T = T(NaN) # Sound speed at the sheath entrance [m s^-1]
    SOL_connection_length::T = T(NaN) # Parallel connection length [m]
    SOL_total_Fx::T = T(NaN) # Total flux expansion [-]
end

function Base.show(io::IO, ::MIME"text/plain", solbox::SOLBox)
    println(io, "Geometry:")
    if solbox.single_null
        println(io, "  Single null")
    else
        println(io, "  Double null")
    end
    println(io, "  Connection length (m): ", solbox.SOL_connection_length)
    println(io, "  Total flux expansion: ", solbox.SOL_total_Fx)

    println(io, "\nIons:")
    println(io, "  [INPUT] - Coefficient of ion conductivity: ", solbox.κ0_i)
    println(io, "  [INPUT] - Ion particle recycling coefficient: ", solbox.recycling_coeff_i)
    println(io, "  [INPUT] - Target ion temperature (eV): ", solbox.Ti_t)
    println(io, "  [INPUT] - Ion mass (amu): ", solbox.mass_ion)
    println(io, "  Upstream parallel ion heat flux (GW m^-2): ", solbox.qpar_i * 1.0e-09)
    println(io, "  Upstream ion particle flux (ptcles m^-2 s^-1): ", solbox.Γ_i)
    println(io, "  Target ion density (ptcls m^-3): ", solbox.ni_t)
    println(io, "  Upstream ion temperature (eV): ", solbox.Ti_u)
    println(io, "  Upstream ion density (ptcls m^-3): ", solbox.ni_u)

    println(io, "\nElectrons:")
    println(io, "  [INPUT] - Coefficient of electron conductivity: ", solbox.κ0_e)
    println(io, "  [INPUT] - Electron particle recycling coefficient: ", solbox.recycling_coeff_e)
    println(io, "  [INPUT] - Target electron temperature (eV): ", solbox.Te_t)
    println(io, "  Upstream parallel electron heat flux (GW m^-2): ", solbox.qpar_e * 1.0e-09)
    println(io, "  Upstream electron particle flux (ptcles m^-2 s^-1): ", solbox.Γ_e)
    println(io, "  Target electron density (ptcls m^-3): ", solbox.ne_t)
    println(io, "  Upstream electron temperature (eV): ", solbox.Te_u)
    println(io, "  Upstream electron density (ptcls m^-3): ", solbox.ne_u)

    println(io, "\nOther quantities:")
    println(io, "  Sound speed (km s^-1): ", solbox.cst * 1.0e-03)
    println(io, "  [INPUT] - frac_cond: ", solbox.frac_cond)
    println(io, "  [INPUT] - frac_mom: ", solbox.frac_mom)
    println(io, "  [INPUT] - λq (mm): ", solbox.λq * 1E3)
    return nothing
end

mutable struct ActorSOLBox{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorSOLBox{P}}
    solbox::SOLBox{D}
end

function ActorSOLBox(dd::IMAS.dd{D}, par::FUSEparameters__ActorSOLBox{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorSOLBox)
    par = OverrideParameters(par; kw...)
    return ActorSOLBox(dd, par, SOLBox{D}())
end

"""
    ActorSOLBox(dd::IMAS.dd, act::ParametersAllActors; kw...)

0D box model for the scrape-off layer developed by X. Zhang et al. https://doi.org/10.1016/j.nme.2022.101354 .
"""
function ActorSOLBox(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorSOLBox(dd, act.ActorSOLBox; kw...)
    step(actor)
    finalize(actor)
    return actor
end

function _step(actor::ActorSOLBox{D,P}) where {D<:Real,P<:Real}
    # Retrieve variables defined in the struct ActorSOLBox
    dd = actor.dd
    par = actor.par

    eqt = dd.equilibrium.time_slice[]
    cp1d = dd.core_profiles.profiles_1d[]

    solbox = actor.solbox
    solbox.κ0_i = par.κ0_i
    solbox.recycling_coeff_i = par.recycling_coeff_i
    solbox.Ti_t = par.Ti_t
    solbox.mass_ion = IMAS.A_effective(cp1d)
    solbox.κ0_e = par.κ0_e
    solbox.recycling_coeff_e = par.recycling_coeff_e
    solbox.Te_t = par.Te_t
    solbox.frac_cond = par.frac_cond
    solbox.frac_mom = par.frac_mom
    solbox.λq = par.λq

    # Step 1 - get the parallel heat and particles fluxes for ions and electrons

    # Retrieve the total ion and electron sources (power and particles)
    total_source = IMAS.total_sources(dd.core_sources, cp1d; time0=dd.global_time, fields=[:power_inside, :total_ion_power_inside, :particles_inside])

    # Retrieve Psol and the particles per second across the separatrix, for both ions and electrons
    power_ions = total_source.total_ion_power_inside[end]
    power_electrons = total_source.electrons.power_inside[end]
    Psol = power_ions + power_electrons
    particles_ions_electrons = total_source.electrons.particles_inside[end] # Same for ions and electrons due to quasineutrality

    # Get the R,Z coordinates of the outer midplane (R = R0 + a)
    r_omp = eqt.profiles_1d.r_outboard[end]
    z_omp = eqt.global_quantities.magnetic_axis.z

    # Get the reference R coordinate and toroidal field
    R0, B0 = eqt.global_quantities.vacuum_toroidal_field.r0, eqt.global_quantities.vacuum_toroidal_field.b0

    # Get the magnetic field at outer midplane
    eqt2d = findfirst(:rectangular, eqt.profiles_2d)
    _, _, PSI_interpolant = IMAS.ψ_interpolant(eqt2d)
    Bp = IMAS.Bp(PSI_interpolant, r_omp, z_omp) # r and z component of B
    Bt = abs.(B0 .* R0 ./ r_omp)                # toroidal component of B
    B = sqrt.(Bp .^ 2 + Bt .^ 2)                # total magnetic field B

    if solbox.λq == Inf
        # Take λq as twice λq as predicted by Eich's Bpol scaling
        solbox.λq = 2.0 * IMAS.widthSOL_eich(eqt, Psol)
    end

    # Calculate the parallel area of the separatrix flux tube between midplane and outer target
    A_par = (2.0 * pi * r_omp * solbox.λq * (Bp / B))

    # Calculate the power and particle flux densities for both ions and electrons
    solbox.qpar_e = power_electrons / A_par
    solbox.qpar_i = power_ions / A_par
    solbox.Γ_e = particles_ions_electrons / A_par
    solbox.Γ_i = solbox.Γ_e # quasineutrality

    # Step 2 - get connection length and total flux expansion for the flux tube

    # Trace field lines in the SOL
    SOL = IMAS.sol(eqt, dd.wall; levels=100)

    # We need to decide which flux surface we will take as being representative of the separatrix. We will not use the
    # exact separatrix here as a) the connection length blows up too close to it, and b) sometimes the flux surface tracing
    # can wrap back around the plasma when it shouldn't, so we want to be a little bit into the SOL (by 'offset' number of flux surfaces).
    offset = 3

    # In order to make this assessment, we also need to check how balanced the two xpoints are - if they are quite balanced then even
    # the last flux tube in the 'near SOL' will actually be right next to the separatrix (which we want to avoid).

    # Get the flux of the primary and secondary separatrices
    LCFS_psi = SOL[:lfs][1].psi
    LCFS2_psi = SOL[:lfs_far][1].psi

    # Check the magnetic balance of the two nulls
    if abs((LCFS_psi - LCFS2_psi) / (abs(LCFS_psi) + abs(LCFS2_psi))) > 0.025
        # Single null
        solbox.single_null = true
        separatrix = SOL[:lfs][1+offset]
    else
        # (Disconnected) double null
        solbox.single_null = false
        separatrix = SOL[:lfs_far][1+offset]
    end

    # Parallel connection length from midplane to target (last point)
    solbox.SOL_connection_length = separatrix.s[end]

    # Total flux expansion from midplane to target (last point)
    solbox.SOL_total_Fx = separatrix.total_flux_expansion[end]

    # Step 3 - calculate the sound speed
    solbox.cst = sqrt(IMAS.mks.e * (solbox.Te_t + solbox.Ti_t) / (solbox.mass_ion * IMAS.mks.m_p))

    # Step 4 - calculate the upstream temperature(s)

    # Intermediate variable. Calculate once for reduced computation.
    alpha = solbox.frac_cond * 1.75 * solbox.SOL_connection_length * (1.0 - 0.5 * (solbox.SOL_total_Fx - 1.0))

    # Upstream ion temperature
    solbox.Ti_u = ((solbox.Ti_t)^(7.0 / 2.0) + alpha * (solbox.qpar_i / solbox.κ0_i))^(2.0 / 7.0)

    # Upstream electron temperature
    solbox.Te_u = ((solbox.Te_t)^(7.0 / 2.0) + alpha * (solbox.qpar_e / solbox.κ0_e))^(2.0 / 7.0)

    # Step 5 - calculate the target densities

    # Calculate the particle transmission coefficients

    # Ion recycling coefficient
    γi = 1.0 - solbox.recycling_coeff_i

    # Electron recycling coefficient
    γe = 1.0 - solbox.recycling_coeff_e

    # Intermediate variable. Calculate once for reduced computation.
    beta = (solbox.SOL_total_Fx + 1) / (2.0 * solbox.SOL_total_Fx)

    # Target ion density
    solbox.ni_t = beta * (solbox.Γ_i / (γi * solbox.cst))

    # Target electron density
    solbox.ne_t = beta * (solbox.Γ_e / (γe * solbox.cst))

    # Step 6 - calculate the upstream densities

    # Upstream ion density
    solbox.ni_u = ((2.0 * solbox.ni_t) / (solbox.frac_mom)) * ((solbox.Ti_t + solbox.Te_t) / (solbox.Ti_u + solbox.Te_u))

    # Upstream electron density
    solbox.ne_u = ((2.0 * solbox.ne_t) / (solbox.frac_mom)) * ((solbox.Ti_t + solbox.Te_t) / (solbox.Ti_u + solbox.Te_u))

    if par.verbose
        display(solbox)
    end

    # Plotting
    if par.do_plot
        # Extract the R, Z of the flux surface (separatrix) trace from outer midplane to target
        lower_outer_target = dd.divertors.divertor[1].target[1].tile[1].surface_outline
        sep_r_plot = separatrix.r[separatrix.midplane_index:end]
        sep_z_plot = separatrix.z[separatrix.midplane_index:end]

        if solbox.single_null
            connection_lengths = [lfs.s[end] for lfs in SOL[:lfs]]
            psi_vals = [lfs.psi for lfs in SOL[:lfs]]
        else
            connection_lengths = [lfs_far.s[end] for lfs_far in SOL[:lfs_far]]
            psi_vals = [lfs_far.psi for lfs_far in SOL[:lfs_far]]
        end

        # Plot
        p1 = plot(; layout=1, size=(500, 700))
        plot!(
            p1,
            SOL;
            colorbar_entry=false,
            color=:grays,
            xlabel="R (m)",
            ylabel="Z (m)",
            title="Connection length: " * string(round(solbox.SOL_connection_length; digits=2)) * " m"
        )
        plot!(p1, sep_r_plot, sep_z_plot; linewidth=2, color=:red, label="Flux tube")
        plot!(p1, lower_outer_target.r, lower_outer_target.z; linewidth=2, color=:blue, label="Outer target")
        p2 = plot(; layout=1, size=(1000, 700))
        plot!(p2, psi_vals, connection_lengths; legend=false, xlabel="ψ", ylabel="L (m)")
        p = plot(p1, p2; layout=(1, 2))
        display(p)
    end

    return actor
end

function _finalize(actor::ActorSOLBox)
    dd = actor.dd

    solbox = actor.solbox

    ep1d = dd.edge_profiles.profiles_1d[]

    # We will populate the edge_profiles IDS - profiles_1d, and, as this model is
    # for just a single flux tube, only 1 point in these 1D radial profiles will be
    # populated

    @assert ep1d.grid.rho_pol_norm[1] == 1.0 "ep1d.grid.rho_pol_norm = $(ep1d.grid.rho_pol_norm)"
    ep1d.electrons.density[1] = solbox.ne_u
    ep1d.electrons.temperature[1] = solbox.Te_u

    # For now, I am only populating the electron density, as this is what will act as a
    # boundary condition for when FUSE solves core transport. In reality, FUSE will actually
    # use the pedestal top density, which can be calculated as something like 4 x the separatrix density.

    return actor
end
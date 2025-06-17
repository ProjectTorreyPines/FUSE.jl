#= ===================== =#
#       ActorSOLBox       #
#= ===================== =#

# Part of act structure (act.ActorSOLBox here)
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
    mass_ion::Entry{T} = Entry{T}("kg", "Ion mass in multiples of amu"; default=2.5)
    recycling_coeff_i::Entry{T} = Entry{T}("-", "Ion particle recycling coefficient"; default=0.99)
    recycling_coeff_e::Entry{T} = Entry{T}("-", "Electron particle recycling coefficient"; default=0.99)
    λ_mm::Entry{T} = Entry{T}("mm", "Width of the flux tube"; default=Inf)
    do_debug::Entry{Bool} = Entry{Bool}("-", "Flag for debugging"; default=false)
    do_plot::Entry{Bool} = act_common_parameters(; do_plot=false)
end

mutable struct ActorSOLBox{D,P} <: SingleAbstractActor{D,P}
    dd::IMAS.dd{D}
    par::OverrideParameters{P,FUSEparameters__ActorSOLBox{P}}
    qpar_i::Union{Nothing,Float64} # Upstream parallel ion heat flux [W m^-2]
    qpar_e::Union{Nothing,Float64} # Upstream parallel electron heat flux [W m^-2]
    Γ_i::Union{Nothing,Float64} # Upstream ion particle flux [ptcles m^-2 s^-1]
    Γ_e::Union{Nothing,Float64}  # Upstream electron particle flux [ptcles m^-2 s^-1]
    Te_u::Union{Nothing,Float64} # Electron temperature at the upstream (separatrix) [eV]
    Ti_u::Union{Nothing,Float64} # Ion temperature at the upstream (separatrix) [eV]
    ne_t::Union{Nothing,Float64} # Electron density at the target [pcles m^-2 s^-1]
    ni_t::Union{Nothing,Float64} # Ion density at the target [pcles m^-2 s^-1]
    ne_u::Union{Nothing,Float64} # Electron density at the upstream (separatrix) [pcles m^-2 s^-1]
    ni_u::Union{Nothing,Float64} # Ion density at the upstream (separatrix) [pcles m^-2 s^-1]
    cst::Union{Nothing,Float64} # Sound speed at the sheath entrance [m s^-1]
    SOL_connection_length::Union{Nothing,Float64} # Parallel connection length [m]
    SOL_total_Fx::Union{Nothing,Float64} # Total flux expansion [-]
end

function ActorSOLBox(dd::IMAS.dd{D}, par::FUSEparameters__ActorSOLBox{P}; kw...) where {D<:Real,P<:Real}
    logging_actor_init(ActorSOLBox)
    par = OverrideParameters(par; kw...)
    return ActorSOLBox(dd, par, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing) # should be the same number of inputs as struct above
end

"""
    ActorSOLBox(dd::IMAS.dd, act::ParametersAllActors; kw...)

    0D box model for the scrape-off layer developed by X. Zhang et al. https://doi.org/10.1016/j.nme.2022.101354 .
"""
# FUSE.ActorSOLBox(dd,act) will trigger: Step -> Fianlize 
function ActorSOLBox(dd::IMAS.dd, act::ParametersAllActors; kw...)
    actor = ActorSOLBox(dd, act.ActorSOLBox; kw...)
    step(actor)
    finalize(actor)
    return actor
end

# Step function
function _step(actor::ActorSOLBox{D,P}) where {D<:Real,P<:Real}

    # Retrieve variables defined in the struct ActorSOLBox
    dd = actor.dd
    par = actor.par

    # Step 1 - get the parallel heat and particles fluxes for ions and electrons

    # Retrieve the total ion and electron sources (power and particles)
    total_source = IMAS.total_sources(dd.core_sources, dd.core_profiles.profiles_1d[]; time0=0.0)

    # Retrieve Psol and the particles per second across the separatrix, for both ions and electrons
    power_ions = total_source.total_ion_power_inside[end]
    power_electrons = total_source.electrons.power_inside[end]
    particles_ions_electrons = total_source.electrons.particles_inside[end] # Same for ions and electrons due to quasineutrality

    # Get the R,Z coordinates of the outer midplane (R = R0 + a)
    eqt = dd.equilibrium.time_slice[]
    r_omp = eqt.profiles_1d.r_outboard[end]
    z_omp = eqt.global_quantities.magnetic_axis.z

    # Get the reference R coordinate and toroidal field
    R0, B0 = eqt.global_quantities.vacuum_toroidal_field.r0, eqt.global_quantities.vacuum_toroidal_field.b0

    eqt2d = findfirst(:rectangular, eqt.profiles_2d)

    # Get the magnetic field at outer midplane
    _, _, PSI_interpolant = IMAS.ψ_interpolant(eqt2d)
    Bp = IMAS.Bp(PSI_interpolant, r_omp, z_omp) # r and z component of B
    Bt = abs.(B0 .* R0 ./ r_omp)                # toroidal component of B
    B = sqrt.(Bp .^ 2 + Bt .^ 2)                # total magnetic field B

    if par.λ_mm == Inf

        # Take λ as twice λq as predicted by Eich's Bpol scaling
        par.λ_mm = 2.0 * 0.63 * (Bp^-1.19)

    end

    # Convert flux tube width mm -> m
    λ = par.λ_mm * 1.0e-03

    # Calculate the parallel area of the separatrix flux tube between midplane and outer target
    A_par = (2.0 * pi * r_omp * λ * (Bp / B))

    # Calculate the power and particle flux densities for both ions and electrons
    actor.qpar_e = power_electrons / A_par
    actor.qpar_i = power_ions / A_par
    actor.Γ_e = particles_ions_electrons / A_par
    actor.Γ_i = actor.Γ_e # quasineutrality

    # Define physical constants
    amu = 1.660538921 * 1.0e-27 # Atomic mass unit [kg]
    e_charge = 1.60e-19 # Magnitude of the electron charge [C]

    # Step 2 - get connection length and total flux expansion for the flux tube

    # Trace field lines in the SOL
    SOL = IMAS.sol(dd; levels=100)

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
    if abs((LCFS_psi - LCFS2_psi) / LCFS_psi) > 0.025

        # Single null
        single_null = true
        separatrix = SOL[:lfs][1+offset]

    else

        # (Disconnected) double null
        single_null = false
        separatrix = SOL[:lfs_far][1+offset]

    end

    # Parallel connection length from midplane to target (last point)
    actor.SOL_connection_length = separatrix.s[end]

    # Total flux expansion from midplane to target (last point)
    actor.SOL_total_Fx = separatrix.total_flux_expansion[end]

    # Step 3 - calculate the sound speed
    actor.cst = sqrt(e_charge * (par.Te_t + par.Ti_t) / (par.mass_ion * amu))

    # Step 4 - calculate the upstream temperature(s)

    # Intermediate variable. Calculate once for reduced computation.
    alpha = par.frac_cond * 1.75 * actor.SOL_connection_length * (1.0 - 0.5 * (actor.SOL_total_Fx - 1.0))

    # Upstream ion temperature
    actor.Ti_u = ((par.Ti_t)^(7.0 / 2.0) + alpha * (actor.qpar_i / par.κ0_i))^(2.0 / 7.0)

    # Upstream electron temperature
    actor.Te_u = ((par.Te_t)^(7.0 / 2.0) + alpha * (actor.qpar_e / par.κ0_e))^(2.0 / 7.0)

    # Step 5 - calculate the target densities

    # Calculate the particle transmission coefficients

    # Ion recycling coefficient
    γi = 1.0 - par.recycling_coeff_i

    # Electron recycling coefficient
    γe = 1.0 - par.recycling_coeff_e

    # Intermediate variable. Calculate once for reduced computation.
    beta = (actor.SOL_total_Fx + 1) / (2.0 * actor.SOL_total_Fx)

    # Target ion density
    actor.ni_t = beta * (actor.Γ_i / (γi * actor.cst))

    # Target electron density
    actor.ne_t = beta * (actor.Γ_e / (γe * actor.cst))

    # Step 6 - calculate the upstream densities

    # Upstream ion density
    actor.ni_u = ((2.0 * actor.ni_t) / (par.frac_mom)) * ((par.Ti_t + par.Te_t) / (actor.Ti_u + actor.Te_u))

    # Upstream electron density
    actor.ne_u = ((2.0 * actor.ne_t) / (par.frac_mom)) * ((par.Ti_t + par.Te_t) / (actor.Ti_u + actor.Te_u))

    if par.do_debug

        println("\nDebugging")

        println("\nGeometry")

        if single_null
            println("Single null")
        else
            println("Double null")
        end

        println("Connection length (m): ", actor.SOL_connection_length)
        println("Total flux expansion: ", actor.SOL_total_Fx)

        println("\nIons")
        println("[INPUT] - Coefficient of ion conductivity: ", par.κ0_i)
        println("[INPUT] - Ion particle recycling coefficient: ", par.recycling_coeff_i)
        println("Upstream parallel ion heat flux (GW m^-2): ", actor.qpar_i * 1.0e-09)
        println("Upstream ion particle flux (ptcles m^-2 s^-1): ", actor.Γ_i)
        println("[INPUT] - Target ion temperature (eV): ", par.Ti_t)
        println("Target ion density (ptcls m^-3): ", actor.ni_t)
        println("Upstream ion temperature (eV): ", actor.Ti_u)
        println("Upstream ion density (ptcls m^-3): ", actor.ni_u)
        println("[INPUT] - Ion mass (amu): ", par.mass_ion)

        println("\nElectrons")
        println("[INPUT] - Coefficient of electron conductivity: ", par.κ0_e)
        println("[INPUT] - Electron particle recycling coefficient: ", par.recycling_coeff_e)
        println("Upstream parallel electron heat flux (GW m^-2): ", actor.qpar_e * 1.0e-09)
        println("Upstream electron particle flux (ptcles m^-2 s^-1): ", actor.Γ_e)
        println("[INPUT] - Target electron temperature (eV): ", par.Te_t)
        println("Target electron density (ptcls m^-3): ", actor.ne_t)
        println("Upstream electron temperature (eV): ", actor.Te_u)
        println("Upstream electron density (ptcls m^-3): ", actor.ne_u)

        println("\nOther quantities")
        println("Sound speed (km s^-1): ", actor.cst * 1.0e-03)
        println("[INPUT] - frac_cond: ", par.frac_cond)
        println("[INPUT] - frac_mom: ", par.frac_mom)
        println("λ (mm): ", par.λ_mm)

    end

    # Plotting
    if par.do_plot

        # Extract the R, Z of the flux surface (separatrix) trace from outer midplane to target
        sep_r = separatrix.r
        sep_z = separatrix.z
        sep_con_len = separatrix.s
        N = length(sep_con_len)
        index_0 = separatrix.midplane_index
        lower_outer_target = dd.divertors.divertor[1].target[1].tile[1].surface_outline
        indeces = [i for i in range(; start=1, step=1, length=N) if i >= index_0]
        sep_r_plot = [sep_r[i] for i in range(; start=1, step=1, length=N) if i >= index_0]
        sep_z_plot = [sep_z[i] for i in range(; start=1, step=1, length=N) if i >= index_0]

        # Get an array of midplane to outer target connection lengths for a set of
        # flux tubes in the SOL
        SOL = IMAS.sol(dd; levels=100)

        if single_null

            Ntubes = length(SOL[:lfs])
            connection_lengths = [SOL[:lfs][i].s[end] for i in range(; start=1, step=1, length=Ntubes)]
            psi_vals = [SOL[:lfs][i].psi for i in range(; start=1, step=1, length=Ntubes)]

        else

            Ntubes = length(SOL[:lfs_far])
            connection_lengths = [SOL[:lfs_far][i].s[end] for i in range(; start=1, step=1, length=Ntubes)]
            psi_vals = [SOL[:lfs_far][i].psi for i in range(; start=1, step=1, length=Ntubes)]

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
            title="Connection length: " * string(round(actor.SOL_connection_length; digits=2)) * " m"
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

# Finalize function
function _finalize(actor::ActorSOLBox)
    # Finalization after computation is done in step
    # (Probably populate dd; for details of dd: https://fuse.help/dev/dd.html)

    # Retrieve variables defined in the struct ActorSOLBox
    dd = actor.dd
    par = actor.par

    # The below is still a work in progress

    # We will populate the edge_profiles IDS - profiles_1d, and, as this model is
    # for just a single flux tube, only 1 point in these 1D radial profiles will be
    # populated

    dd.edge_profiles.profiles_1d[].electrons.density[1] = actor.ne_u
    dd.edge_profiles.profiles_1d[].electrons.temperature[1] = actor.Te_u

    # For now, I am only populating the electron density, as this is what will act as a
    # boundary condition for when FUSE solves core transport. In reality, FUSE will actually
    # use the pedestal top density, which can be calculated as something like 4 x the separatrix density.

    return actor
end